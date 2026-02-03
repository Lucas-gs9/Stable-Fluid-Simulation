#include <algorithm>
#include <array>
#include <fstream>
#include <random>
#include <execution>
#include<iostream>

#include "Labs/Final_proj/CaseStableFluids.h"
#include "Labs/Common/ImGuiHelper.h"

using namespace VCX::Labs::Final_proj;

static constexpr auto c_Size = std::pair(1024U, 1024U);

int FluidState::neighbor(int x, int dx, int max) const {
    if (wrapBorder) {
        x += dx;
        return (x % max + max) % max;
    } else {
        return std::clamp(x + dx, 0, max - 1);
    }
}

float FluidState::neighbor(float x, int dx, int max) const {
    if (wrapBorder) {
        return std::fmod(x + float(dx) + float(max), float(max));
    } else {
        return std::clamp(x + float(dx), 0.f, float(max - 1));
    }
}

void FluidState::Resize(Size2D const& s, int v)
{
    size          = s;
    if (v == 0) {
        std::size_t n = static_cast<std::size_t>(size.first) * size.second;

        uCur.resize(n);
        uNext.resize(n);
        vCur.resize(n);
        vNext.resize(n);
        densityrCur.resize(n);
        densityrNext.resize(n);
        densitygCur.resize(n);
        densitygNext.resize(n);
        densitybCur.resize(n);
        densitybNext.resize(n);
        div.resize(n);
        p.resize(n);
    }
    InitbgFlow();
}

void FluidState::SwapBuffers()
{
    uCur.swap(uNext);
    vCur.swap(vNext);
    densityrCur.swap(densityrNext);
    densitygCur.swap(densitygNext);
    densitybCur.swap(densitybNext);
}

void FluidState::Transport(std::vector<float>& S1, const std::vector<float>& S0, const std::vector<float>& U, const std::vector<float>& V, float dt,int i,int j) 
{
    const int Nx = size.first;
    const int Ny = size.second;
    auto      idx = [=](int x, int y) { return y * Nx + x; };
    float     u   = U[idx(i, j)];
    float     v   = V[idx(i, j)];
    float     x0  = float(i) - u * dt;
    float     y0  = float(j) - v * dt;
    x0            = neighbor(x0, 0, Nx);
    y0            = neighbor(y0, 0, Ny);
    // bilinear interpolation
    int   ix      = int(x0);
    int   iy      = int(y0);
    float fx      = x0 - ix;
    float fy      = y0 - iy;
    int   ix1     = (ix + 1) % Nx;
    int   iy1     = (iy + 1) % Ny;
    float s00     = S0[idx(ix, iy)];
    float s01     = S0[idx(ix, iy1)];
    float s10     = S0[idx(ix1, iy)];
    float s11     = S0[idx(ix1, iy1)];
    S1[idx(i, j)] = (s00 * (1 - fx) + s10 * fx) * (1 - fy) + (s01 * (1 - fx) + s11 * fx) * fy;
}

void FluidState::Diffuse(std::vector<float>& S1, const std::vector<float>& S0, const std::vector<float>& oS0, float k, float dt, int i, int j) 
{
    const int Nx = size.first;
    const int Ny  = size.second;
    auto      idx = [=](int x, int y) { return y * Nx + x; };

    int   iL  = neighbor(i, -1, Nx);
    int   iR  = neighbor(i, 1, Nx);
    int   jM  = neighbor(j, -1, Ny);
    int   jP  = neighbor(j, 1, Ny);
    float sum = S0[idx(iL, j)] + S0[idx(iR, j)] + S0[idx(i, jM)] + S0[idx(i, jP)];

    S1[idx(i, j)] = (oS0[idx(i, j)] + dt * k * sum) / (1.f + 4 * dt * k);
}

void FluidState::Project(std::vector<float>& U1, std::vector<float>& V1, const std::vector<float>& U0, const std::vector<float>& V0)
{
    const int Nx  = size.first;
    const int Ny  = size.second;
    const int N   = Nx * Ny;
    auto      idx = [=](int x, int y) { return y * Nx + x; };
    //calculate div
    std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(N), [&](int t) {
        int   i    = t % Nx;
        int   j    = t / Nx;
        int   iL   = neighbor(i, -1, Nx);
        int   iR   = neighbor(i, 1, Nx);
        int   jM   = neighbor(j, -1, Ny);
        int   jP   = neighbor(j, 1, Ny);
        float dudx = U0[idx(iR, j)] - U0[idx(iL, j)];
        float dvdy = V0[idx(i, jP)] - V0[idx(i, jM)];
        div[t]     = 0.5f * (dudx + dvdy);
    });
    //Poisson
    std::fill(p.begin(), p.end(), 0.0f);
    std::vector<float> pOld = p;
    const int          iter = 15;
    for (int k = 0; k < iter; ++k) {
        std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(N), [&](int t) {
            int   i   = t % Nx;
            int   j   = t / Nx;
            int   iL  = neighbor(i, -1, Nx);
            int   iR  = neighbor(i, 1, Nx);
            int   jM  = neighbor(j, -1, Ny);
            int   jP  = neighbor(j, 1, Ny);
            float sum = pOld[idx(iL, j)] + pOld[idx(iR, j)]
                + pOld[idx(i, jM)] + pOld[idx(i, jP)];
            p[t] = (sum - div[t]) * 0.25f;
        });
        std::swap(p, pOld);
    }
    if (iter & 1) std::copy(pOld.begin(), pOld.end(), p.begin());

    std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(N), [&](int t) {
        int   i    = t % Nx;
        int   j    = t / Nx;
        int   iL   = neighbor(i, -1, Nx);
        int   iR   = neighbor(i, 1, Nx);
        int   jM   = neighbor(j, -1, Ny);
        int   jP   = neighbor(j, 1, Ny);
        float dpdx = 0.5f * (p[idx(iR, j)] - p[idx(iL, j)]);
        float dpdy = 0.5f * (p[idx(i, jP)] - p[idx(i, jM)]);
        U1[t]      = U0[t] - dpdx;
        V1[t]      = V0[t] - dpdy;
    });
}

void FluidState::Dissipate(std::vector<float>& S1, const std::vector<float>& S0, float a, float dt,int i)
{
    const float  factor = 1.0f - std::min(a * dt, 1.0f);
    S1[i]              = S0[i] * factor; 
}

CaseStableFluids::CaseStableFluids() :
    _texture(Engine::GL::SamplerOptions {
        .MinFilter = Engine::GL::FilterMode::Linear,
        .MagFilter = Engine::GL::FilterMode::Nearest }),
    _front(Common::CreatePureImageRGB(512, 512, { 0, 0, 0 })){
    ResetFluid();
}

void CaseStableFluids::OnSetupPropsUI() {
    ImGui::Checkbox("Zoom Tooltip", &_enableZoom);
    ImGui::ColorEdit3("Tint", &_fluid._tint[0], ImGuiColorEditFlags_NoInputs);
    ImGui::Checkbox("Infinite Domain", &_fluid.wrapBorder);
    ImGui::SliderFloat("Amount", &_fluid.denAmount, 0.0f, 2.0f, "%.2f");
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("0 = pure velocity disturb");
    ImGui::SliderInt("Brush Radius", &_brushRadius, 1, 7);
    if (ImGui::Button("Reset")) 
    {
        _running = false;
        _task.Reset();
        ResetFluid();
        _needStep = true; 
        _fluid.ReloadbgFlow();
    }
    static const char * items[] = { "128", "256", "512", "1024" };
    static int          gridId  = 3;
    if (ImGui::Combo("Grid", &gridId, items, IM_ARRAYSIZE(items))) {
        int newN = 128 << gridId;
        if (newN != static_cast<int>(_fluid.size.first*2)) {
            _running = false;
            _task.Reset();
            _fluid.Resize({ static_cast<uint32_t>(newN / 2), static_cast<uint32_t>(newN / 2) });
            ResetFluid();
            _needStep = true; 
            _fluid.ReloadbgFlow();
        }
    }

    ImGui::SliderFloat("DT", &_fluid.dt, 0.01f, 0.5f, "%.3f");
    ImGui::SliderFloat("Diffusion", &_fluid.diff, 0.f, 0.01f, "%.5f");
    ImGui::SliderFloat("Viscosity", &_fluid.visc, 0.f, 0.01f, "%.5f");
    ImGui::SliderFloat("Density Fade", &_fluid.densityDissipation, 0.0f, 0.1f, "%.3f");
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("0 = keep forever, 0.1 = vanish quickly");
    ImGui::SliderFloat("Brightness", &_fluid.brightness, 0.1f, 2.f, "%.05f");
    const char * bgNames[] = { "None", "Horizontal", "Vertical", "Vortex" };
    int          bg        = static_cast<int>(_fluid._bgFlowType);
    if (ImGui::Combo("Background Flow", &bg, bgNames, IM_ARRAYSIZE(bgNames))) {
        _running = false;
        _task.Reset();
        ResetFluid();
        _needStep          = true; 
        _fluid._bgFlowType = static_cast<BgFlowType>(bg);
        _fluid.ReloadbgFlow();
    }
}

VCX::Labs::Common::CaseRenderResult CaseStableFluids::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
    auto const [width, height] = _fluid.size;

    if (_front.GetSizeX() != width * 2 || _front.GetSizeY() != height * 2) {
        _front = Common::CreatePureImageRGB(width * 2, height * 2, { 0, 0, 0 });
    }

    float dt = ImGui::GetIO().DeltaTime;
    _accPhys += dt;

    if ((_accPhys >= _physDt || _needStep) && ! _running) {
        _accPhys  = 0.f;
        _needStep = false;
        _task.Emplace([this]() {
            StepFluid(); 
            Common::ImageRGB img(_fluid.size.first * 2, _fluid.size.second * 2);
            RGBFromDensity(img);
            return img;
        });
        _running = true;
    }
    if (_running && _task.HasValue()) 
    {
        _front   = _task.Value();
        _running = false;
        _needStep = true;
    }
    _texture.Update(_task.ValueOr(_front));
    
    return Common::CaseRenderResult {
        .Fixed     = true,
        .Image     = _texture,
        .ImageSize = {_fluid.size.first * 2, _fluid.size.second * 2},
    };
}

void CaseStableFluids::OnProcessInput(ImVec2 const& pos) {
    if (_enableZoom) {
        auto window  = ImGui::GetCurrentWindow();
        bool hovered = false, held = false;
        ImGui::ButtonBehavior(window->Rect(), window->GetID("##io"), &hovered, &held);
        if (! held && hovered)
            Common::ImGuiHelper::ZoomTooltip(_texture, { _fluid.size.first * 2, _fluid.size.second * 2 }, pos);
    }

    ImGuiIO &  io            = ImGui::GetIO();
    glm::ivec2        curPos = { pos.x, pos.y };
    const auto [simW, simH]   = _fluid.size; 
    const glm::ivec2 dispSize = { int(2 * simW), int(2 * simH) };
    bool       inside = curPos.x >= 0 && curPos.y >= 0 && curPos.x < dispSize.x && curPos.y < dispSize.y;

    static glm::ivec2 lastPos(-1);
    if (io.MouseDown[ImGuiMouseButton_Left] && inside) {
        if (lastPos.x >= 0) {
            glm::vec2 lastSim = glm::vec2(lastPos) * 0.5f;
            glm::vec2 curSim  = glm::vec2(curPos) * 0.5f;

            int Nx = _fluid.size.first, Ny = _fluid.size.second;
            int steps          = 1 + int(glm::length(curSim - lastSim) / 0.4f); 
            steps              = glm::clamp(steps, 3, 8 + Nx / 32);
            glm::vec2 deltaSim = glm::vec2(io.MouseDelta.x, io.MouseDelta.y) * 0.5f;

            for (int step = 1; step <= steps; ++step) {
                glm::vec2 p = glm::mix(lastSim, curSim, float(step) / float(steps));
                if (p.x < 0 || p.x >= Nx || p.y < 0 || p.y >= Ny) continue;
                glm::vec2 vel = deltaSim * 20.f;

                int R  = _brushRadius;
                int cx = int(p.x), cy = int(p.y);
                for (int dy = -R; dy <= R; ++dy)
                    for (int dx = -R; dx <= R; ++dx) {
                        int x = cx + dx;
                        int y = cy + dy;
                        if (_fluid.wrapBorder) {
                            x = (x + Nx) % Nx;
                            y = (y + Ny) % Ny;
                        } else {
                            if (x < 0 || x >= Nx || y < 0 || y >= Ny) continue;
                        }
                        float r2 = float(dx * dx + dy * dy) / (R * R);
                        if (r2 > 1.f) continue;
                        float fall = 1.f - pow(r2, 0.25);

                        AddDensity(x, y, _fluid.denAmount * fall);
                        AddForce(x, y, vel.x * fall, vel.y * fall);
                    }
            }
        }
        lastPos = curPos;
    } else {
        lastPos = glm::ivec2(-1);
    }
}

void CaseStableFluids::ResetFluid() {
    if (_fluid.size.first == 0 || _fluid.size.second == 0) {
        _fluid.Resize({ 512, 512 }, 0);
    }
    std::size_t n = static_cast<std::size_t>(_fluid.size.first) * _fluid.size.second;
    std::fill(_fluid.uCur.begin(), _fluid.uCur.end(), 0.f);
    std::fill(_fluid.uNext.begin(), _fluid.uNext.end(), 0.f);
    std::fill(_fluid.vCur.begin(), _fluid.vCur.end(), 0.f);
    std::fill(_fluid.vNext.begin(), _fluid.vNext.end(), 0.f);
    std::fill(_fluid.densityrCur.begin(), _fluid.densityrCur.end(), 0.f);
    std::fill(_fluid.densityrNext.begin(), _fluid.densityrNext.end(), 0.f);
    std::fill(_fluid.densitygCur.begin(), _fluid.densitygCur.end(), 0.f);
    std::fill(_fluid.densitygNext.begin(), _fluid.densitygNext.end(), 0.f);
    std::fill(_fluid.densitybCur.begin(), _fluid.densitybCur.end(), 0.f);
    std::fill(_fluid.densitybNext.begin(), _fluid.densitybNext.end(), 0.f);
    std::fill(_fluid.div.begin(), _fluid.div.end(), 0.f);
    std::fill(_fluid.p.begin(), _fluid.p.end(), 0.f);
    _front = Common::CreatePureImageRGB(_fluid.size.first * 2, _fluid.size.second * 2, glm::u8vec3 { 0, 0, 0 });
    _texture.Update(_front);
    _needStep = true;
}

void CaseStableFluids::StepFluid() {
    const int Nx   = _fluid.size.first;
    const int Ny   = _fluid.size.second;
    float     dt   = _fluid.dt;
    float     visc = _fluid.visc;
    float     diff = _fluid.diff;
    int       iter = 15;
    int       N    = Nx * Ny;

    std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(N), [&](int idx) {
        int i = idx % Nx;
        int j = idx / Nx;
        _fluid.Transport(_fluid.densityrNext, _fluid.densityrCur, _fluid.uCur, _fluid.vCur, dt, i, j);
        _fluid.Transport(_fluid.densitygNext, _fluid.densitygCur, _fluid.uCur, _fluid.vCur, dt, i, j);
        _fluid.Transport(_fluid.densitybNext, _fluid.densitybCur, _fluid.uCur, _fluid.vCur, dt, i, j);
        _fluid.Transport(_fluid.uNext, _fluid.uCur, _fluid.uCur, _fluid.vCur, dt, i, j);
        _fluid.Transport(_fluid.vNext, _fluid.vCur, _fluid.uCur, _fluid.vCur, dt, i, j);
    }); // v: cur->next, den:cur->next
    
    std::vector<float> odenr = _fluid.densityrNext;
    std::vector<float> odeng = _fluid.densitygNext;
    std::vector<float> odenb = _fluid.densitybNext;

    std::vector<float> ou   = _fluid.uNext;
    std::vector<float> ov   = _fluid.vNext;
    for (int iter_k=0;iter_k<iter;++iter_k) {
        std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(N), [&](int idx) {
            int i = idx % Nx;
            int j = idx / Nx;
            _fluid.Diffuse(_fluid.densityrCur, _fluid.densityrNext, odenr, diff, dt, i, j);
            _fluid.Diffuse(_fluid.densitygCur, _fluid.densitygNext, odeng, diff, dt, i, j);
            _fluid.Diffuse(_fluid.densitybCur, _fluid.densitybNext, odenb, diff, dt, i, j);
            _fluid.Diffuse(_fluid.uCur, _fluid.uNext, ou, visc, dt, i, j);
            _fluid.Diffuse(_fluid.vCur, _fluid.vNext, ov, visc, dt, i, j);
        });
        std::swap(_fluid.densityrCur, _fluid.densityrNext);
        std::swap(_fluid.densitygCur, _fluid.densitygNext);
        std::swap(_fluid.densitybCur, _fluid.densitybNext);
        std::swap(_fluid.uCur, _fluid.uNext);
        std::swap(_fluid.vCur, _fluid.vNext);
    }//v:next->cur, den:next->cur

    _fluid.Project(_fluid.uNext, _fluid.vNext, _fluid.uCur, _fluid.vCur);//v:cur->next

     std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(N), [&](int i) {
        _fluid.Dissipate(_fluid.densityrCur, _fluid.densityrCur, _fluid.densityDissipation, dt, i);
        _fluid.Dissipate(_fluid.densitygCur, _fluid.densitygCur, _fluid.densityDissipation, dt, i);
        _fluid.Dissipate(_fluid.densitybCur, _fluid.densitybCur, _fluid.densityDissipation, dt, i);
        _fluid.Dissipate(_fluid.uCur, _fluid.uNext, 0.01f, dt, i);
        _fluid.Dissipate(_fluid.vCur, _fluid.vNext, 0.01f, dt, i);
        _fluid.uCur[i] += _fluid._bgFlow[_fluid._bgFlowType][i][0] * std::min(0.01f * dt, 1.f);
        _fluid.vCur[i] += _fluid._bgFlow[_fluid._bgFlowType][i][1] * std::min(0.01f * dt, 1.f);
    }); // v:next->cur, den: cur->cur
}

void CaseStableFluids::AddDensity(int x, int y,float amount) {
    auto [W, H] = _fluid.size;
    if (W == 0 || H == 0) return;

    int    ix  = std::clamp(x, 0, static_cast<int>(W) - 1);
    int    iy  = std::clamp(y, 0, static_cast<int>(H) - 1);
    size_t idx = iy * W + ix;

    _fluid.densityrCur[idx] += _fluid._tint.r * amount;
    _fluid.densitygCur[idx] += _fluid._tint.g * amount;
    _fluid.densitybCur[idx] += _fluid._tint.b * amount;
}

void CaseStableFluids::AddForce(int x, int y, float fx, float fy) {
    auto [W, H] = _fluid.size;
    if (W == 0 || H == 0) return;

    int    ix  = std::clamp(x, 0, static_cast<int>(W) - 1);
    int    iy  = std::clamp(y, 0, static_cast<int>(H) - 1);
    size_t idx = iy * W + ix;

    _fluid.uCur[idx] += fx;
    _fluid.vCur[idx] += fy;
}

void CaseStableFluids::RGBFromDensity(Common::ImageRGB & img) {
    const auto [simW, simH] = _fluid.size;
    const int dispW         = simW * 2;
    const int dispH         = simH * 2;

    const float displayCap = 4.0f;
    const float brightness = _fluid.brightness;

    std::vector<glm::vec3> res(simW * simH);
    std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(simW * simH), [&](int t) {
        glm::vec3 d = { _fluid.densityrCur[t], _fluid.densitygCur[t], _fluid.densitybCur[t] };
        glm::vec3 v = glm::clamp(d / displayCap, 0.0f, 1.0f);
        res[t]     = { v.r * std::pow(brightness, 0.75), v.g * std::pow(brightness, 0.75), v.b * std::pow(brightness, 0.75) };
    });

    auto sample = [&](float x, float y) -> glm::vec3 {
        x        = std::clamp(x, 0.0f, float(simW - 1));
        y        = std::clamp(y, 0.0f, float(simH - 1));
        int   ix = int(x), iy = int(y);
        float fx = x - ix, fy = y - iy;
        int   ix1 = std::min(ix + 1, int(simW - 1));
        int   iy1 = std::min(iy + 1, int(simH - 1));
        glm::vec3 r00 = res[iy * simW + ix];
        glm::vec3 r01 = res[iy1 * simW + ix];
        glm::vec3 r10 = res[iy * simW + ix1];
        glm::vec3 r11 = res[iy1 * simW + ix1];
        return (r00 * (1 - fx) + r10 * fx) * (1 - fy) + (r01 * (1 - fx) + r11 * fx) * fy;
    };

    std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(dispW * dispH), [&](int t) {
        int di = t % dispW;
        int dj = t / dispW;

        float x        = (di + 0.5f) / 2.0f - 0.5f;
        float y        = (dj + 0.5f) / 2.0f - 0.5f;
        glm::vec3 g        = sample(x, y);
        img.At(di, dj) = g;
    });
}

void FluidState::InitbgFlow()
{
    const int Nx = size.first;
    const int Ny = size.second;
    const int N  = Nx * Ny;
    for (auto & table : _bgFlow) table.resize(N);
    std::fill(_bgFlow[BG_NONE].begin(), _bgFlow[BG_NONE].end(), glm::vec2(0.f));

    std::fill(_bgFlow[BG_HORZ].begin(), _bgFlow[BG_HORZ].end(), glm::vec2(10.f, 0.f));

    std::fill(_bgFlow[BG_VERT].begin(), _bgFlow[BG_VERT].end(), glm::vec2(0.f, -10.f));

    glm::vec2 c(Nx * 0.5f, Ny * 0.5f);
    std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(N), [&](int idx) {
        int       i = idx % Nx, j = idx / Nx;
        glm::vec2 p(i + 0.5f, j + 0.5f);
        glm::vec2 dir = p - c;
        float     r   = glm::length(dir);
        glm::vec2 v(0.f);
        if (r > 0.1f) {
            glm::vec2 tan(-dir.y, dir.x);
            v = glm::normalize(tan) * std::exp(-r * 0.005f) * 8.f;
        }
        _bgFlow[BG_VORTEX][idx] = v;
    });
}

void FluidState::ReloadbgFlow()
{
    const int Nx = size.first;
    const int Ny = size.second;
    const int N  = Nx * Ny;
    std::for_each(std::execution::par_unseq, counting_iterator<int>(0), counting_iterator<int>(N), [&](int idx) {
        uCur[idx] += _bgFlow[_bgFlowType][idx][0];
        vCur[idx] += _bgFlow[_bgFlowType][idx][1];
    });
}