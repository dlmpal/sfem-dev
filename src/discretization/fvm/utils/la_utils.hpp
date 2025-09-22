#pragma once

#include "../fv_function.hpp"
#include "../../../la/petsc/sfem_petsc.hpp"

namespace sfem::fvm
{
    /// @brief Create a vector for a finite volume function
    la::Vector create_vec(const FVFunction &phi);

    /// @brief Create a matrix for a finite volume function
    la::SparseMatrix create_mat(const FVFunction &phi);
}

#ifdef SFEM_HAS_PETSC

namespace sfem::fvm::petsc
{
    /// @brief Create a PETSc vector for a finite volume function
    la::petsc::PetscVec create_vec(const FVFunction &phi);

    /// @brief Create a PETSc matrix for a finite volume function
    la::petsc::PetscMat create_mat(const FVFunction &phi);
}

#endif // SFEM_HAS_PETSC