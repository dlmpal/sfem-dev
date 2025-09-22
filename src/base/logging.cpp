#include "logging.hpp"
#include "application.hpp"
#include "../parallel/mpi.hpp"
#include <format>

namespace sfem
{
    //=============================================================================
    static std::string log_level_str(LogLevel level)
    {
        switch (level)
        {
        case LogLevel::debug:
            return std::string("debug");
        case LogLevel::info:
            return std::string("info");
        case LogLevel::warning:
            return std::string("warning");
        default:
            return std::string("error");
        }
    }
    //=============================================================================
    void log_msg(const std::string &msg,
                 bool root_only,
                 LogLevel level,
                 std::source_location location)
    {
        // Format the message
        std::string formatted;
        formatted += std::format("[{}]-[{}]-[{}]: ", Application::instance().name(),
                                 mpi::rank(),
                                 log_level_str(level));
        formatted += msg;

        // Warnings and error also includ file and line information
        if (level >= LogLevel::warning)
        {
            formatted += std::format("\t at file {}, line {}\n", location.file_name(), location.line());
        }

        if ((root_only and mpi::rank() == mpi::root()) or !root_only)
        {
            Application::instance().log_message(formatted, level);
        }

        // For errors, abort the application
        if (level == LogLevel::error)
        {
            mpi::abort();
        }
    }
}