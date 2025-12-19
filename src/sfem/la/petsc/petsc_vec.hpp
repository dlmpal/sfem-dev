#pragma once

#ifdef SFEM_HAS_PETSC

#include <sfem/base/config.hpp>
#include <petscvec.h>
#include <vector>
#include <span>

namespace sfem::la::petsc
{
    /// @brief Thin wrapper around PETSc's Vec
    class PetscVec
    {
    public:
        /// @brief Initialize from an existing Vec
        /// @param x Existing Vec
        /// @param inc_ref_count Whether to increase the reference count for the PETSc object
        PetscVec(Vec x, bool inc_ref_count);

        // Avoid copies
        PetscVec(const PetscVec &) = delete;
        PetscVec &operator=(const PetscVec &) = delete;

        // Move constructor and assignment
        PetscVec(PetscVec &&);
        PetscVec &operator=(PetscVec &&);

        // Destructor
        ~PetscVec();

        /// @brief Get the local size
        int size_local() const;

        /// @brief Get the global size
        int size_global() const;

        /// @brief Get the underlying Vec
        Vec vec() const;

        /// @brief Copy the vector
        PetscVec copy() const;

        /// @brief Set all vector values
        void set_all(real_t value);

        /// @brief Set values into the vector
        /// @param idxs Indices
        /// @param values Values
        /// @param mode Whether to insert or add the values
        void set_values(std::span<const int> idxs,
                        std::span<const real_t> values,
                        InsertMode mode = ADD_VALUES);

        /// @brief Assemble the vector
        void assemble();

        /// @brief Get the (local) values
        /// @note The values for ghost indices are also included
        std::vector<real_t> get_values() const;

    private:
        /// @brief Underlying PETSc Vec
        Vec vec_;
    };
}

#endif // SFEM_HAS_PETSC