#include "Assets/bundled.h"
#include "Labs/Final_proj/App.h"
#include<iostream>

int main() {
    using namespace VCX;
    return Engine::RunApp<Labs::Final_proj::App>(Engine::AppContextOptions {
        .Title      = "VCX Final Project",
        .WindowSize = { 1024, 768 },
        .FontSize   = 16,

        .IconFileNames = Assets::DefaultIcons,
        .FontFileNames = Assets::DefaultFonts,
    });
}
