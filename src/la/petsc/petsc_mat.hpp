#pragma once

#ifdef SFEM_HAS_PETSC

#include "petsc.hpp"

namespace sfem::la::petsc
{
    /// @brief Thin wrapper around the PETSc Mat
    class PetscMat
    {
    public:
        /// @brief Create a PetscMat from an existing PETSc Mat
        /// @param A Existing PETSc Mat
        /// @param inc_ref_count Whether to increase the ref count for A
        PetscMat(Mat mat, bool inc_ref_count);

        // Copy constructor (deleted)
        PetscMat(const PetscMat &) = delete;

        // Copy assignment operator (deleted)
        PetscMat &operator=(const PetscMat &) = delete;

        /// @brief Move constructor
        PetscMat(PetscMat &&);

        /// @brief Move assignment
        PetscMat &operator=(PetscMat &&);

        /// @brief Destructor
        ~PetscMat();

        /// @brief Get the local size (no. rows and cols)
        std::array<int, 2> size_local() const;

        /// @brief Get the global size (no. rows and cols)
        std::array<int, 2> size_global() const;

        /// @brief Get the underlying PETSc Mat
        Mat mat() const;

        /// @brief Reset the matrix's memory
        /// @note Call before re-assembling
        void reset();

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