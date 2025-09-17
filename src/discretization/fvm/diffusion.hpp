#pragma once

#include "fv_bc.hpp"
#include "../coefficient.hpp"
#include "../../la/native/setval_utils.hpp"

namespace sfem::fvm
{
    void diffusion(const FVFunction &phi,
                   const FVFunction &grad,
                   const FVBC &bc,
                   const Coefficient &coeff,
                   la::MatSet lhs, la::VecSet rhs);
}