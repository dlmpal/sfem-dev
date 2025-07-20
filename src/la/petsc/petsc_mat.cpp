#ifdef SFEM_HAS_PETSC

#include "petsc_mat.hpp"
#include "../../base/error.hpp"

namespace sfem::la::petsc
{
    //=============================================================================
    PetscMat::PetscMat(Mat mat, bool inc_ref_count)
        : mat_(mat)
    {
        if (!mat_)
        {
            SFEM_ERROR("Invalid PETSc Mat passed to PetscMat constructor");
        }

        if (inc_ref_count)
        {
            PetscObjectReference((PetscObject)mat_);
        }
    }
    //=============================================================================
    PetscMat::PetscMat(PetscMat &&other)
    {
        mat_ = other.mat_;
        other.mat_ = nullptr;
    }
    //=============================================================================
    PetscMat &PetscMat::operator=(PetscMat &&other)
    {
        if (this != &other)
        {
            if (mat_)
            {
                MatDestroy(&mat_);
            }
            mat_ = other.mat_;
            other.mat_ = nullptr;
        }

        return *this;
    }
    //=============================================================================
    PetscMat::~PetscMat()
    {
        if (mat_)
        {
            MatDestroy(&mat_);
        }
    }
    //=============================================================================
    std::array<int, 2> PetscMat::size_local() const
    {
        int n_rows, n_cols;
        MatGetLocalSize(mat_, &n_rows, &n_cols);
        return {n_rows, n_cols};
    }
    //=============================================================================
    std::array<int, 2> PetscMat::size_global() const
    {
        int n_rows, n_cols;
        MatGetSize(mat_, &n_rows, &n_cols);
        return {n_rows, n_cols};
    }
    //=============================================================================
    Mat PetscMat::mat() const
    {
        return mat_;
    }
    //=============================================================================
    void PetscMat::reset()
    {
        MatResetPreallocation(mat_);
    }
    //=============================================================================
    void PetscMat::set_values(std::span<const int> row_idxs,
                              std::span<const int> col_idxs,
                              std::span<const real_t> values)
    {
        // Get the block size and check that sizes match
        int block_size;
        MatGetBlockSize(mat_, &block_size);
        SFEM_CHECK_SIZES(row_idxs.size() * col_idxs.size() * block_size * block_size, values.size());

        if (block_size == 1)
        {
            MatSetValuesLocal(mat_,
                              static_cast<int>(row_idxs.size()),
                              row_idxs.data(),
                              static_cast<int>(col_idxs.size()),
                              col_idxs.data(),
                              values.data(),
                              ADD_VALUES);
        }
        else
        {
            MatSetValuesBlockedLocal(mat_,
                                     static_cast<int>(row_idxs.size()),
                                     row_idxs.data(),
                                     static_cast<int>(col_idxs.size()),
                                     col_idxs.data(),
                                     values.data(),
                                     ADD_VALUES);
        }
    }
    //=============================================================================
    void PetscMat::assemble()
    {
        MatAssemblyBegin(mat_, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat_, MAT_FINAL_ASSEMBLY);
    }
}

#endif // SFEM_HAS_PETSC