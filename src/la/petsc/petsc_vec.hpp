#pragma once

#ifdef SFEM_HAS_PETSC

#include "petsc.hpp"

namespace sfem::la::petsc
{
    /// @brief Thin wrapper around the PETSc Vec
    class PetscVec
    {
    public:
        /// @brief Create a PetscVec from an existing PETSc Vec
        /// @param x Existing PETSc Vec
        /// @param inc_ref_count Whether to increase the reference count for x
        PetscVec(Vec x, bool inc_ref_count);

        // Copy constructor (deleted)
        PetscVec(const PetscVec &) = delete;

        // Copy assignment operator (deleted)
        PetscVec &operator=(const PetscVec &) = delete;

        /// @brief Move constructor
        PetscVec(PetscVec &&);

        /// @brief Move assignment
        PetscVec &operator=(PetscVec &&);

        /// @brief Destructor
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
        /// @param insert Whether to insert the values, overriding previous
        void set_values(std::span<const int> idxs,
                        std::span<const real_t> values,
                        bool insert = false);

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