#pragma once

#include "function.hpp"

namespace sfem
{
    class Coefficient
    {
    public:
        virtual real_t &operator()(int idx, int comp) = 0;
        virtual real_t operator()(int idx, int comp) const = 0;
        virtual real_t operator()(const std::array<real_t, 3> &pt, real_t time) const = 0;
        virtual real_t &operator()(const std::array<real_t, 3> &pt, real_t time) = 0;
    };

    class ConstantCoefficient : public Coefficient
    {
    public:
        ConstantCoefficient(real_t value);

        real_t &operator()(int idx, int comp) override;
        real_t operator()(int idx, int comp) const override;
        real_t operator()(const std::array<real_t, 3> &pt, real_t time) const override;
        real_t &operator()(const std::array<real_t, 3> &pt, real_t time) override;

    private:
        /// @brief Constant value
        real_t value_;
    };
}