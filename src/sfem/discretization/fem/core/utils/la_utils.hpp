#pragma once

#include <sfem/discretization/fem/core/fe_field.hpp>
#include <sfem/la/backend.hpp>
#include <sfem/la/petsc/sfem_petsc.hpp>

namespace sfem::fem
{
    /// @brief Create a vector for a finite element field
    la::Vector create_vec(const FEField &phi);

    /// @brief Create a matrix for a finite element field
    la::SparseMatrix create_mat(const FEField &phi);

    /// @brief Create a linear system for a finite element field
    std::shared_ptr<la::LinearSystem> create_axb(const FEField &phi,
                                                 la::SolverType solver_type = la::SolverType::gmres,
                                                 la::SolverOptions solver_options = {},
                                                 la::Backend backend = la::Backend::native);
}

#ifdef SFEM_HAS_PETSC

namespace sfem::fem::petsc
{
    /// @brief Create a PETSc vector for a finite element field
    la::petsc::PetscVec create_vec(const FEField &phi);

    /// @brief Create a PETSc matrix for a finite element field
    la::petsc::PetscMat create_mat(const FEField &phi);
}

#endif // SFEM_HAS_PETSC