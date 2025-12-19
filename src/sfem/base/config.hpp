#pragma once

namespace sfem
{
    #ifdef SFEM_USE_SINGLE_PRECISION
        using real_t = float;
    #elif defined SFEM_USE_DOUBLE_PRECISION
        using real_t = double;
    #endif
}