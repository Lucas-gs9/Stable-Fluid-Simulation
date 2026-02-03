#pragma once

#include <vector>

#include "Engine/app.h"
#include "Labs/Final_proj/CaseStableFluids.h"
#include "Labs/Common/UI.h"

namespace VCX::Labs::Final_proj {
    class App : public Engine::IApp {
    private:
        Common::UI              _ui;
        std::size_t             _caseId = 0;
        CaseStableFluids            _caseStableFluids;

        std::vector<std::reference_wrapper<Common::ICase>> _cases = {
            _caseStableFluids
        };

    public:
        App();

        void OnFrame() override;
    };
}
#pragma once
