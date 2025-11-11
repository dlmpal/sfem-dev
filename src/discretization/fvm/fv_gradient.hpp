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

    /// @brief Create the gradient field for a given finite volume field
    /// @note Does not actually compute the gradient
    /// @param phi Finite volume field
    /// @return Gradient field
    [[nodiscard]] std::shared_ptr<FVField> create_gradient(const FVField &phi);

    /// @brief Compute the gradient of a finite volume field
    /// @param phi Finite volume field
    /// @param bc Boundary condition
    /// @param grad Field gradient
    /// @param method Gradient evaluation method
    void gradient(const FVField &phi,
                  const FVBC &bc,
                  FVField &grad,
                  GradientMethod method = GradientMethod::green_gauss);

}