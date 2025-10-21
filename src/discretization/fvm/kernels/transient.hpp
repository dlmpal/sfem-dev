#pragma once

#include "../fv_field.hpp"
#include "../../coefficient.hpp"
#include "../../../la/native/setval_utils.hpp"

namespace sfem::fvm
{
    /// @brief Assemble the LHS matrix and RHS vector for a first-order temporal derivative
    /// acting on a scalar finite volume field
    /// @param phi Finite volume field
    /// @param coeff Cell coefficient
    /// @param dt Timestep size
    /// @param lhs LHS matrix
    /// @param rhs RHS vector
    void ddt(const FVField &phi, const Coefficient &coeff,
             real_t dt, la::MatSet lhs, la::VecSet rhs);
}