#include "Labs/Final_proj/App.h"

namespace VCX::Labs::Final_proj {
    App::App():
        _ui(Labs::Common::UIOptions {}) {
    }

    void App::OnFrame() {
        _ui.Setup(_cases, _caseId);
    }
}