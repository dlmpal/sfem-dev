#pragma once

#include <sfem/base/application.hpp>

namespace sfem
{
    /// @brief Initializes SFEM and returns the singleton application instance.
    /// Should be called at the start of every SFEM program
    Application &initialize(int argc, char *argv[],
                            const std::string &app_name = "sfem",
                            const std::filesystem::path &log_filename = {});
}