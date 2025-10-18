#pragma once

#include "../fv_bc.hpp"
#include "../../coefficient.hpp"
#include "../../../la/native/setval_utils.hpp"

namespace sfem::fvm
{
    /// @brief Assemble the LHS matrix and RHS vector for the diffusion operator
    /// acting on a scalar finite volume field
    /// @param phi Finite volume field
    /// @param grad Field gradient
    /// @param bc Boundary condition
    /// @param dt Timestep size
    /// @param cell_coeff Cell coefficient
    /// @param facet_coeff Facet coefficient
    /// @param lhs LHS matrix
    /// @param rhs RHS vector
    void diffusion(const FVField &phi,
                   const FVField &grad,
                   const FVBC &bc, real_t dt,
                   const Coefficient &cell_coeff,
                   const Coefficient &facet_coeff,
                   la::MatSet lhs, la::VecSet rhs);
}