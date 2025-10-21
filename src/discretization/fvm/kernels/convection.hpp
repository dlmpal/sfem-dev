#pragma once

#include "../fv_bc.hpp"
#include "../../coefficient.hpp"
#include "../../../la/native/setval_utils.hpp"

namespace sfem::fvm
{
    enum class DifferencingScheme
    {
        upwind,
        linear,
        linear_upwind,
        blended
    };

    /// @brief Assemble the LHS matrix and RHS vector for the convenction operator
    /// acting on a scalar finite volume field
    /// @param phi Finite volume field
    /// @param grad Field gradient
    /// @param bc Boundary condition
    /// @param vel Convecting velocity field
    /// @param lhs LHS matrix
    /// @param rhs RHS vector
    /// @param scheme Differencing scheme
    /// @param implicit Whether to treat this term implicitly (add contributions to LHS),
    /// or explicitly (add contributions to RHS)
    void convection(const FVField &phi,
                    const FVField &grad,
                    const FVBC &bc,
                    const Coefficient &vel,
                    la::MatSet lhs, la::VecSet rhs,
                    DifferencingScheme scheme = DifferencingScheme::upwind,
                    bool implicit = true);
}