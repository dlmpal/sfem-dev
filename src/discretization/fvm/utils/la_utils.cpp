#include "la_utils.hpp"

namespace sfem::fvm::petsc
{
    //=============================================================================
    la::petsc::PetscVec create_vec(const FVFunction &phi)
    {
        const auto fv_space = phi.space();
        return la::petsc::create_vec(*fv_space->index_map(), phi.n_comp());
    }
    //=============================================================================
    la::petsc::PetscMat create_mat(const FVFunction &phi)
    {
        const auto fv_space = phi.space();
        return la::petsc::create_mat(*fv_space->connectivity(),
                                     *fv_space->index_map(),
                                     *fv_space->index_map(),
                                     phi.n_comp());
    }
}