#include "la_utils.hpp"
#include "../../../la/native/sparse_matrix.hpp"

namespace sfem::fem
{
    //=============================================================================
    la::Vector create_vec(const FEFunction &phi)
    {
        const auto fe_space = phi.space();
        return la::Vector(fe_space->index_map(), phi.n_comp());
    }
    //=============================================================================
    la::SparseMatrix create_mat(const FEFunction &phi)
    {
        const auto fe_space = phi.space();
        return la::SparseMatrix(fe_space->connectivity()[1],
                                fe_space->index_map(),
                                fe_space->index_map(),
                                phi.n_comp());
    }
}

#ifdef SFEM_HAS_PETSC

namespace sfem::fem::petsc
{
    //=============================================================================
    la::petsc::PetscVec create_vec(const FEFunction &phi)
    {
        const auto fe_space = phi.space();
        return la::petsc::create_vec(*fe_space->index_map(), phi.n_comp());
    }
    //=============================================================================
    la::petsc::PetscMat create_mat(const FEFunction &phi)
    {
        const auto fe_space = phi.space();
        return la::petsc::create_mat(*fe_space->connectivity()[1],
                                     *fe_space->index_map(),
                                     *fe_space->index_map(),
                                     phi.n_comp());
    }
    //=============================================================================
    void solve(la::petsc::PetscMat &A,
               la::petsc::PetscVec &b,
               la::petsc::PetscVec &x,
               const DirichletBC &bc)
    {
        auto [vars, values] = bc.get_dofs_values();
        la::petsc::eliminate_rows_cols(vars, values,
                                       A, b, x);
        la::petsc::solve(A, b, x);
    }
}

#endif // SFEM_HAS_PETSC