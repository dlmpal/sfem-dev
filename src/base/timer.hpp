#pragma once

#include "logging.hpp"
#include <chrono>

namespace sfem
{
    /// @brief Time the execution of a function
    class Timer
    {
    public:
        Timer(std::source_location location = std::source_location::current());
        ~Timer();

    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> start_;
        std::chrono::time_point<std::chrono::high_resolution_clock> stop_;
        std::source_location location_;
    };
}