#pragma once

#include "logging.hpp"
#include <format>

#define SFEM_ERROR(msg)                                   \
    {                                                     \
        sfem::log_msg(msg, false, sfem::LogLevel::error); \
    }

#define SFEM_CHECK_FILE_OPEN(file, filename)                                       \
    {                                                                              \
        if (file.is_open() == false)                                               \
        {                                                                          \
            auto msg = std::format("Could not open file {}\n", filename.string()); \
            SFEM_ERROR(msg);                                                       \
        }                                                                          \
    }

#define SFEM_CHECK_SIZES(good, bad)                                                    \
    {                                                                                  \
        if (static_cast<int>(good) != static_cast<int>(bad))                           \
        {                                                                              \
            auto msg = std::format("Got size {} while expected size {}\n", bad, good); \
            SFEM_ERROR(msg);                                                           \
        }                                                                              \
    }

#define SFEM_CHECK_INDEX(idx, range)                                                  \
    {                                                                                 \
        if (idx < 0 or idx >= range)                                                  \
        {                                                                             \
            auto msg = std::format("Index {} is out of range [0, {})\n", idx, range); \
            SFEM_ERROR(msg);                                                          \
        }                                                                             \
    }
