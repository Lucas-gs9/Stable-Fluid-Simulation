#pragma once

#include "Engine/Async.hpp"
#include "Labs/Common/ICase.h"
#include "Labs/Common/ImageRGB.h"

namespace VCX::Labs::Final_proj {
    using Size2D = std::pair<std::uint32_t, std::uint32_t>;

    enum BgFlowType : std::size_t {
        BG_NONE = 0,
        BG_HORZ,
        BG_VERT,
        BG_VORTEX,
        BG_COUNT
    };

    template<class T>
    struct counting_iterator {
        T value {};
        counting_iterator() = default;
        explicit counting_iterator(T v):
            value(v) {}
        T                   operator*() const { return value; }
        counting_iterator & operator++() {
            ++value;
            return *this;
        }
        bool operator!=(counting_iterator rhs) const { return value != rhs.value; }

        counting_iterator & operator+=(std::ptrdiff_t n) {
            value += n;
            return *this;
        }
        counting_iterator operator+(std::ptrdiff_t n) const {
            auto tmp = *this;
            tmp += n;
            return tmp;
        }
        std::ptrdiff_t operator-(counting_iterator rhs) const { return value - rhs.value; }
        bool           operator<(counting_iterator rhs) const { return value < rhs.value; }
    };

    template<class T>
    struct std::iterator_traits<counting_iterator<T>> {
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using iterator_category = std::random_access_iterator_tag;
    };

    struct FluidState {
        Size2D size { 0, 0 };
        glm::vec3 _tint { 1.0f, 1.0f, 1.0f };
        std::vector<float> uCur, uNext;
        std::vector<float> vCur, vNext; // velocity
        std::vector<float> densityrCur, densityrNext;
        std::vector<float> densitygCur, densitygNext;
        std::vector<float> densitybCur, densitybNext;
        std::vector<float> div, p;                  // divergence, pressure
        std::array<std::vector<glm::vec2>, BG_COUNT> _bgFlow;

        float dt   = 0.1f; // time
        float diff = 0.0001f;//diffusitivity
        float visc = 0.0001f;//viscosity
        float denAmount              = 1.5f;
        float densityDissipation = 0.02f;
        float brightness             = 1.2f;
        bool  wrapBorder         = true;
        BgFlowType _bgFlowType            = BG_NONE; 

        int  neighbor(int x, int dx, int max) const;
        float neighbor(float x, int dx, int max) const;
        void Resize(Size2D const & s, int v = 1);
        void SwapBuffers();
        void  Transport(std::vector<float> & S1, const std::vector<float> & S0, const std::vector<float> & U, const std::vector<float> & V, float dt, int i, int j);
        void  Diffuse(std::vector<float> & S1, const std::vector<float> & S0, const std::vector<float>&oS0, float k, float dt, int i, int j);
        void Project(std::vector<float> & U1, std::vector<float> & V1, const std::vector<float> & U0, const std::vector<float> & V0);
        void  Dissipate(std::vector<float> & S1, const std::vector<float> & S0, float a, float dt, int i);
        void  InitbgFlow();
        void  ReloadbgFlow();
    };

    class CaseStableFluids : public Common::ICase {
    public:
        CaseStableFluids();
        virtual std::string_view const GetName() override { return "Stable Fluids"; }
        virtual void                   OnSetupPropsUI() override;
        virtual Common::CaseRenderResult OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) override;
        virtual void                     OnProcessInput(ImVec2 const & pos) override;

    private:
        FluidState                  _fluid;
        Engine::GL::UniqueTexture2D _texture;
        Common::ImageRGB            _front;
        Engine::Async<Common::ImageRGB> _task;

        bool _enableZoom = true;
        float       _accPhys    = 0.f;
        const float _physDt     = 0.032f;
        bool        _needStep   = false;
        bool _running    = false;
        int       _brushRadius = 5;

        void ResetFluid();
        void StepFluid(); 
        void AddDensity(int x, int y, float amount);
        void AddForce(int x, int y, float fx, float fy);
        void RGBFromDensity(Common::ImageRGB & img);
    };
}
