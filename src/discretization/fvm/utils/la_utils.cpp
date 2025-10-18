#include "la_utils.hpp"
#include "../../../la/native/sparse_matrix.hpp"

namespace sfem::fvm
{
    //=============================================================================
    la::Vector create_vec(const FVField &phi)
    {
        const auto fv_space = phi.space();
        return la::Vector(fv_space->index_map(), phi.n_comp());
    }
    //=============================================================================
    la::SparseMatrix create_mat(const FVField &phi)
    {
        const auto fv_space = phi.space();
        return la::SparseMatrix(fv_space->connectivity(),
                                fv_space->index_map(),
                                fv_space->index_map(),
                                phi.n_comp());
    }
}

#ifdef SFEM_HAS_PETSC

namespace sfem::fvm::petsc
{
    //=============================================================================
    la::petsc::PetscVec create_vec(const FVField &phi)
    {
        const auto fv_space = phi.space();
        return la::petsc::create_vec(*fv_space->index_map(), phi.n_comp());
    }
    //=============================================================================
    la::petsc::PetscMat create_mat(const FVField &phi)
    {
        const auto fv_space = phi.space();
        return la::petsc::create_mat(*fv_space->connectivity(),
                                     *fv_space->index_map(),
                                     *fv_space->index_map(),
                                     phi.n_comp());
    }
}

#endif // SFEM_HAS_PETSC