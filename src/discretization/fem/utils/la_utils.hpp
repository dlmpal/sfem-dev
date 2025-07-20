#pragma once

#include "assembly.hpp"
#include "../dirichlet_bc.hpp"
#include "../../../la/petsc/sfem_petsc.hpp"

#ifdef SFEM_HAS_PETSC

namespace sfem::fem::petsc
{
    /// @brief Create a PETSc matrix for a finite element function
    la::petsc::PetscMat create_mat(const FEFunction &phi);

    /// @brief Create a PETSc vector for a finite element function
    la::petsc::PetscVec create_vec(const FEFunction &phi);

    /// @brief Form and solve the linear system Ax=b.
    /// Before the system is solved, the Dirichlet boundary condition is applied
    void solve(la::petsc::PetscMat &A,
               la::petsc::PetscVec &b,
               la::petsc::PetscVec &x,
               const DirichletBC &bc);
}

#endif // SFEM_HAS_PETSC