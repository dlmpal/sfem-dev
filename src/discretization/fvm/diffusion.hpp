#pragma once

#include "fv_bc.hpp"
#include "../coefficient.hpp"
#include "../../la/native/setval_utils.hpp"

namespace sfem::fvm
{
    /// @brief Assemble the LHS matrix and RHS vector for a diffusion operator
    /// acting on a scalar finite-volume function
    /// @param phi Finite-volume function
    /// @param grad Function gradient
    /// @param bc Boundary condition
    /// @param coeff Diffusion coefficient
    /// @param lhs LHS matrix
    /// @param rhs RHS vector
    void diffusion(const FVFunction &phi,
                   const FVFunction &grad,
                   const FVBC &bc,
                   const Coefficient &coeff,
                   la::MatSet lhs, la::VecSet rhs);
}