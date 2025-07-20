#pragma once

#include "../base/config.hpp"
#include <functional>
#include <span>

namespace sfem::la
{
    using VecSet = std::function<void(std::span<const int>,
                                      std::span<const real_t>)>;

    using MatSet = std::function<void(std::span<const int>,
                                      std::span<const int>,
                                      std::span<const real_t>)>;
}