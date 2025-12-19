#include "timer.hpp"
#include <format>

namespace sfem
{
    //=============================================================================
    Timer::Timer(LogLevel level, std::source_location location)
        : level_(level),
          location_(location)
    {
        start_ = std::chrono::high_resolution_clock::now();
    }
    //=============================================================================
    Timer::~Timer()
    {
        stop_ = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_ - start_);
        std::string msg = std::format("{} completed in {} milliseconds\n", location_.function_name(), duration.count());
        log_msg(msg, true, level_);
    }
}