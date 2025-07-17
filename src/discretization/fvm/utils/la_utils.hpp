#pragma once

#ifdef SFEM_HAS_PETSC

#include "../fv_function.hpp"
#include "../../../la/petsc/sfem_petsc.hpp"

namespace sfem::fvm::petsc
{
    la::petsc::PetscMat create_mat(const FVFunction &phi);

    la::petsc::PetscVec create_vec(const FVFunction &phi);
}

#endif // SFEM_HAS_PETSC