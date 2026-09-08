#pragma once

#include <nlohmann/json.hpp>

#define SFEM_CHECK_JSON_ENTRY(entry, j)                                                     \
    {                                                                                       \
        if (j.contains(entry) == false)                                                     \
        {                                                                                   \
            auto msg = std::format("Missing required entry for {} in input file\n", entry); \
            SFEM_ERROR(msg);                                                                \
        }                                                                                   \
    }