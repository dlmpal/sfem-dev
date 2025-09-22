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
    /// @param msg Message
    /// @param root_only Whether to output the message on the root process only
    /// @param level Log level
    /// @param location Callsite
    void log_msg(const std::string &msg,
                 bool root_only = false,
                 LogLevel level = LogLevel::info,
                 std::source_location location = std::source_location::current());
}