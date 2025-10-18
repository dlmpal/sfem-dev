#pragma once

#include "fv_field.hpp"
#include "../../la/native/setval_utils.hpp"

namespace sfem::fvm
{
    void add_source_term(const fvm::FVField &phi, la::VecSet vecset);
}