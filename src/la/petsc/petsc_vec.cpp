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
            PetscObjectReference((PetscObject)vec_);
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
                VecDestroy(&vec_);
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
            VecDestroy(&vec_);
        }
    }
    //=============================================================================
    int PetscVec::size_local() const
    {
        int n;
        VecGetLocalSize(vec_, &n);
        return n;
    }
    //=============================================================================
    int PetscVec::size_global() const
    {
        int n;
        VecGetSize(vec_, &n);
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
        VecDuplicate(vec_, &vec);
        VecCopy(vec_, vec);
        PetscVec copy(vec, true);
        return copy;
    }
    //=============================================================================
    void PetscVec::set_all(real_t value)
    {
        VecSet(vec_, value);
    }
    //=============================================================================
    void PetscVec::set_values(std::span<const int> idxs,
                              std::span<const real_t> values,
                              InsertMode mode)
    {
        // Get the block size
        int block_size;
        VecGetBlockSize(vec_, &block_size);

        /// @todo Why is this required?
        if (mode == INSERT_VALUES)
        {
            block_size = 1;
        }

        SFEM_CHECK_SIZES(idxs.size() * block_size, values.size());

        if (block_size == 1)
        {
            VecSetValuesLocal(vec_,
                              static_cast<int>(idxs.size()),
                              idxs.data(),
                              values.data(),
                              mode);
        }
        else
        {
            VecSetValuesBlockedLocal(vec_,
                                     static_cast<int>(idxs.size()),
                                     idxs.data(),
                                     values.data(),
                                     mode);
        }
    }
    //=============================================================================
    void PetscVec::assemble()
    {
        VecAssemblyBegin(vec_);
        VecAssemblyEnd(vec_);
    }
    //=============================================================================
    std::vector<real_t> PetscVec::get_values() const
    {
        VecGhostUpdateBegin(vec_,
                            INSERT_VALUES,
                            SCATTER_FORWARD);
        VecGhostUpdateEnd(vec_,
                          INSERT_VALUES,
                          SCATTER_FORWARD);

        Vec localform; ///< Contains ghost padding
        VecGhostGetLocalForm(vec_, &localform);

        int local_size; ///< Local size (+ ghosts)
        VecGetLocalSize(localform, &local_size);

        /// Copy the values from the local form's array
        real_t *localform_array;
        VecGetArray(localform, &localform_array);
        std::vector<real_t> values(local_size);
        std::copy(localform_array,
                  localform_array + local_size,
                  values.data());

        /// Restore array and local form
        VecRestoreArray(localform, &localform_array);
        VecGhostRestoreLocalForm(vec_, &localform);

        return values;
    }
}

#endif // SFEM_HAS_PETSC