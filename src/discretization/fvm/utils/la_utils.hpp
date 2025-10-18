#pragma once

#include "../fv_field.hpp"
#include "../../../la/petsc/sfem_petsc.hpp"

namespace sfem::fvm
{
    /// @brief Create a vector for a finite volume field
    la::Vector create_vec(const FVField &phi);

    /// @brief Create a matrix for a finite volume field
    la::SparseMatrix create_mat(const FVField &phi);
}

#ifdef SFEM_HAS_PETSC

namespace sfem::fvm::petsc
{
    /// @brief Create a PETSc vector for a finite volume field
    la::petsc::PetscVec create_vec(const FVField &phi);

    /// @brief Create a PETSc matrix for a finite volume field
    la::petsc::PetscMat create_mat(const FVField &phi);
}

#endif // SFEM_HAS_PETSC