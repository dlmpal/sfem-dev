#pragma once

#include "fv_bc.hpp"

namespace sfem::fvm
{
    /// @brief Available gradient evaluation methods
    enum GradientMethod
    {
        green_gauss,
        least_squares
    };

    /// @brief Compute the gradient of a finite-volume function
    /// @param phi Finite-volume function
    /// @param bc Boundary condition
    /// @param grad Function gradient
    /// @param method Gradient evaluation method
    void gradient(const FVFunction &phi,
                  const FVBC &bc,
                  FVFunction &grad,
                  GradientMethod method = GradientMethod::green_gauss);

    /// @brief Compute the gradient of a finite-volume function
    /// @param phi Finite-volume function
    /// @param bc Boundary condition
    /// @param method Gradient evaluation method
    /// @return Function gradient
    [[nodiscard]] std::shared_ptr<FVFunction> gradient(const FVFunction &phi,
                                                       const FVBC &bc,
                                                       GradientMethod method = GradientMethod::green_gauss);

}