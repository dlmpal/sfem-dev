#pragma once

#include <sfem/base/application.hpp>

namespace sfem
{
    /// @brief Initializes SFEM and returns the singleton application instance.
    /// Should be called at the start of every SFEM program
    Application &initialize(int argc, char *argv[],
                            bool write_log_file = false,
                            const std::string &app_name = "sfem");
}