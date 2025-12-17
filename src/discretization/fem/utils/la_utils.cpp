#include "la_utils.hpp"
#include "../../../la/native/sparse_matrix.hpp"

namespace sfem::fem
{
    //=============================================================================
    la::Vector create_vec(const FEField &phi)
    {
        const auto fe_space = phi.space();
        return la::Vector(fe_space->index_map(), phi.n_comp());
    }
    //=============================================================================
    la::SparseMatrix create_mat(const FEField &phi)
    {
        const auto fe_space = phi.space();
        return la::SparseMatrix(fe_space->connectivity()[1],
                                fe_space->index_map(),
                                fe_space->index_map(),
                                phi.n_comp());
    }
    //=============================================================================
    std::shared_ptr<la::LinearSystem> create_axb(const FEField &phi,
                                                 la::SolverType solver_type,
                                                 la::SolverOptions solver_options,
                                                 la::Backend backend)
    {
        const auto im = phi.space()->index_map();
        const auto conn = phi.space()->connectivity().back();
        const int n_comp = phi.n_comp();

        switch (backend)
        {
        case la::Backend::native:
            return std::make_shared<la::NativeLinearSystem>(im, conn,
                                                            solver_type,
                                                            solver_options,
                                                            n_comp);
        case la::Backend::petsc:
#ifdef SFEM_HAS_PETSC
            return std::make_shared<la::petsc::PetscLinearSystem>(*im, *conn,
                                                                  solver_type,
                                                                  solver_options,
                                                                  n_comp);
#else
            log_msg("SFEM has not been compiled with PETSc.\n Falling back to native LA backend\n");
            return create_axb(phi, solver_type, solver_options, la::Backend::native);
#endif
        default:
            break;
        }
        return nullptr;
    }
}

#ifdef SFEM_HAS_PETSC

namespace sfem::fem::petsc
{
    //=============================================================================
    la::petsc::PetscVec create_vec(const FEField &phi)
    {
        const auto fe_space = phi.space();
        return la::petsc::create_vec(*fe_space->index_map(), phi.n_comp());
    }
    //=============================================================================
    la::petsc::PetscMat create_mat(const FEField &phi)
    {
        const auto fe_space = phi.space();
        return la::petsc::create_mat(*fe_space->connectivity()[1],
                                     *fe_space->index_map(),
                                     *fe_space->index_map(),
                                     phi.n_comp());
    }
}

#endif // SFEM_HAS_PETSC