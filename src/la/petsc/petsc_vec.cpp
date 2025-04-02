#ifdef SFEM_HAS_PETSC

#include "petsc_vec.hpp"
#include "../../base/error.hpp"

namespace sfem::la::petsc
{
    //=============================================================================
    PetscVec::PetscVec(Vec vec, bool inc_ref_count) : vec_(vec)
    {
        if (!vec_)
        {
            SFEM_ERROR("Invalid PETSc Vec passed to PetscVec constructor");
        }

        if (inc_ref_count)
        {
            SFEM_CHECK_PETSC_ERROR(PetscObjectReference((PetscObject)vec_));
        }
    }
    //=============================================================================
    PetscVec::PetscVec(PetscVec &&other)
    {
        vec_ = other.vec_;
        other.vec_ = nullptr;
    }
    //=============================================================================
    PetscVec &PetscVec::operator=(PetscVec &&other)
    {
        if (this != &other)
        {
            if (vec_)
            {
                SFEM_CHECK_PETSC_ERROR(VecDestroy(&vec_));
            }
            vec_ = other.vec_;
            other.vec_ = nullptr;
        }

        return *this;
    }
    //=============================================================================
    PetscVec::~PetscVec()
    {
        if (vec_)
        {
            SFEM_CHECK_PETSC_ERROR(VecDestroy(&vec_));
        }
    }
    //=============================================================================
    int PetscVec::size_local() const
    {
        int n;
        SFEM_CHECK_PETSC_ERROR(VecGetLocalSize(vec_, &n));
        return n;
    }
    //=============================================================================
    int PetscVec::size_global() const
    {
        int n;
        SFEM_CHECK_PETSC_ERROR(VecGetSize(vec_, &n));
        return n;
    }
    //=============================================================================
    Vec PetscVec::vec() const
    {
        return vec_;
    }
    //=============================================================================
    PetscVec PetscVec::copy() const
    {
        Vec vec;
        SFEM_CHECK_PETSC_ERROR(VecDuplicate(vec_, &vec));
        SFEM_CHECK_PETSC_ERROR(VecCopy(vec_, vec));
        PetscVec copy(vec, true);
        return copy;
    }
    //=============================================================================
    void PetscVec::set_all(real_t value)
    {
        SFEM_CHECK_PETSC_ERROR(VecSet(vec_, value));
    }
    //=============================================================================
    void PetscVec::set_values(std::span<const int> idxs,
                              std::span<const real_t> values,
                              bool insert)
    {
        // Get the block size
        int block_size;
        SFEM_CHECK_PETSC_ERROR(VecGetBlockSize(vec_, &block_size));

        // Whether to insert (thereby overriding) or add the values
        InsertMode mode = ADD_VALUES;
        if (insert)
        {
            mode = INSERT_VALUES;
            block_size = 1;
        }

        SFEM_CHECK_SIZES(idxs.size() * block_size, values.size());

        int error_code;
        if (block_size == 1)
        {
            error_code = VecSetValues(vec_,
                                      static_cast<int>(idxs.size()),
                                      idxs.data(),
                                      values.data(),
                                      mode);
        }
        else
        {
            error_code = VecSetValuesBlocked(vec_,
                                             static_cast<int>(idxs.size()),
                                             idxs.data(),
                                             values.data(),
                                             mode);
        }
        SFEM_CHECK_PETSC_ERROR(error_code);
    }
    //=============================================================================
    void PetscVec::assemble()
    {
        SFEM_CHECK_PETSC_ERROR(VecAssemblyBegin(vec_));
        SFEM_CHECK_PETSC_ERROR(VecAssemblyEnd(vec_));
    }
    //=============================================================================
    std::vector<real_t> PetscVec::get_values() const
    {
        SFEM_CHECK_PETSC_ERROR(VecGhostUpdateBegin(vec_,
                                                   INSERT_VALUES,
                                                   SCATTER_FORWARD));
        SFEM_CHECK_PETSC_ERROR(VecGhostUpdateEnd(vec_,
                                                 INSERT_VALUES,
                                                 SCATTER_FORWARD));

        Vec localform; ///< Contains ghost padding
        SFEM_CHECK_PETSC_ERROR(VecGhostGetLocalForm(vec_, &localform));

        int local_size; ///< Local size (+ ghosts)
        SFEM_CHECK_PETSC_ERROR(VecGetLocalSize(localform, &local_size));

        /// Copy the values from the local form's array
        real_t *localform_array;
        SFEM_CHECK_PETSC_ERROR(VecGetArray(localform, &localform_array));
        std::vector<real_t> values(local_size);
        std::copy(localform_array,
                  localform_array + local_size,
                  values.data());

        /// Restore array and local form
        SFEM_CHECK_PETSC_ERROR(VecRestoreArray(localform, &localform_array));
        SFEM_CHECK_PETSC_ERROR(VecGhostRestoreLocalForm(vec_, &localform));

        return values;
    }
}

#endif // SFEM_HAS_PETSC