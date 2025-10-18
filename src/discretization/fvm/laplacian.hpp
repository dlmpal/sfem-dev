#pragma once

#include "fv_bc.hpp"
#include "../coefficient.hpp"
#include "../../la/native/setval_utils.hpp"

namespace sfem::fvm
{
    /// @brief Assemble the LHS matrix and RHS vector for the Laplacian operator
    /// acting on a scalar finite volume field
    /// @param phi Finite volume field
    /// @param grad Field gradient
    /// @param bc Boundary condition
    /// @param coeff Facet coefficient
    /// @param lhs LHS matrix
    /// @param rhs RHS vector
    void laplacian(const FVField &phi,
                   const FVField &grad,
                   const FVBC &bc,
                   const Coefficient &coeff,
                   la::MatSet lhs, la::VecSet rhs);
}