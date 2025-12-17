#pragma once

namespace sfem::fvm
{
    /// @brief Available gradient evaluation methods
    enum class GradientMethod
    {
        none,
        green_gauss,
        least_squares
    };

    // Forward declaration;
    class FVField;

    void green_gauss_gradient(FVField phi);

    void least_squares_gradient(FVField phi);
}