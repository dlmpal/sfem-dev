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

    /// @brief Compute the gradient of a finite volume field
    /// @param phi Finite volume field
    /// @param bc Boundary condition
    /// @param grad Field gradient
    /// @param method Gradient evaluation method
    void gradient(const FVField &phi,
                  const FVBC &bc,
                  FVField &grad,
                  GradientMethod method = GradientMethod::green_gauss);

    /// @brief Compute the gradient of a finite volume field
    /// @param phi Finite volume field
    /// @param bc Boundary condition
    /// @param method Gradient evaluation method
    /// @return Field gradient
    [[nodiscard]] std::shared_ptr<FVField> gradient(const FVField &phi,
                                                    const FVBC &bc,
                                                    GradientMethod method = GradientMethod::green_gauss);

}