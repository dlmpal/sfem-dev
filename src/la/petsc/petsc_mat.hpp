#pragma once

#ifdef SFEM_HAS_PETSC

#include "../../base/config.hpp"
#include <petscmat.h>
#include <span>

namespace sfem::la::petsc
{
    /// @brief Thin wrapper around PETSc's Mat
    class PetscMat
    {
    public:
        /// @brief Initialize from an existing Mat
        /// @param A Existing Mat
        /// @param inc_ref_count Whether to increase the ref count for the PETSc object
        PetscMat(Mat mat, bool inc_ref_count);

        // Avoid copies
        PetscMat(const PetscMat &) = delete;
        PetscMat &operator=(const PetscMat &) = delete;

        // Move constructor and assignment
        PetscMat(PetscMat &&);
        PetscMat &operator=(PetscMat &&);

        // Destructor
        ~PetscMat();

        /// @brief Get the local size (no. rows and cols)
        std::array<int, 2> size_local() const;

        /// @brief Get the global size (no. rows and cols)
        std::array<int, 2> size_global() const;

        /// @brief Get the underlying Mat
        Mat mat() const;

        /// @brief Zero all matrix entries
        /// @note Call before re-assembling
        void zero_entries();

        /// @brief Set values into the matrix
        /// @param row_idxs Row indices
        /// @param col_idxs Column indices
        /// @param values Values to be inserted
        void set_values(std::span<const int> row_idxs,
                        std::span<const int> col_idxs,
                        std::span<const real_t> values);

        /// @brief Assemble the matrix
        void assemble();

    private:
        /// @brief Underlying PETSc Mat
        Mat mat_;
    };
}

#endif // SFEM_HAS_PETSC