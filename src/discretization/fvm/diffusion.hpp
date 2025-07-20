#pragma once

#include "fv_function.hpp"
#include "../coefficient.hpp"
#include "../../la/utils.hpp"

namespace sfem::fvm
{
    void isotropic_diffusion(const FVFunction &phi,
                             const Coefficient &coeff,
                             const std::string &region,
                             la::MatSet lhs, la::VecSet rhs);

    void dirichlet_bc(const FVFunction &phi,
                      const Coefficient &coeff,
                      const std::string &region,
                      real_t value,
                      la::MatSet lhs, la::VecSet rhs);
}