#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/la/backend.hpp>
#include <sfem/la/native/linear_system.hpp>
#include <sfem/la/petsc/sfem_petsc.hpp>

namespace sfem::fvm
{
    /// @brief Create a vector for a finite volume field
    la::Vector create_vec(const FVField &phi);

    /// @brief Create a matrix for a finite volume field
    la::SparseMatrix create_mat(const FVField &phi);

    /// @brief Create a linear system for a finite volume field
    std::shared_ptr<la::LinearSystem> create_axb(const FVField &phi,
                                                 la::SolverType solver_type = la::SolverType::gmres,
                                                 la::SolverOptions solver_options = {},
                                                 la::Backend backend = la::Backend::native);
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