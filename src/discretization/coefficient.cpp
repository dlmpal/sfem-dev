#include "coefficient.hpp"

namespace sfem
{
    //=============================================================================
    ConstantCoefficient::ConstantCoefficient(real_t value)
        : value_(value)
    {
    }
    //=============================================================================
    real_t &ConstantCoefficient::operator()([[maybe_unused]] int,
                                            [[maybe_unused]] int)
    {
        return value_;
    }
    //=============================================================================
    real_t ConstantCoefficient::operator()([[maybe_unused]] int,
                                           [[maybe_unused]] int) const
    {
        return value_;
    }
    //=============================================================================
    real_t ConstantCoefficient::operator()([[maybe_unused]] const std::array<real_t, 3> &,
                                           [[maybe_unused]] real_t) const
    {
        return value_;
    }
    //=============================================================================
    real_t &ConstantCoefficient::operator()([[maybe_unused]] const std::array<real_t, 3> &,
                                            [[maybe_unused]] real_t)
    {
        return value_;
    }
}