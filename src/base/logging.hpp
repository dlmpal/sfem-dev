#pragma once

#include <string>
#include <source_location>

namespace sfem
{
    enum class LogLevel
    {
        debug = 0,
        info = 1,
        warning = 2,
        error = 3
    };

    /// @brief Log a message
    void log_msg(const std::string &msg,
                 LogLevel level,
                 std::source_location location = std::source_location::current());
}