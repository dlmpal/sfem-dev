#pragma once

#ifdef SFEM_HAS_PETSC

#include "petsc.hpp"

namespace sfem::la::petsc
{
    /// @brief Thin wrapper around the PETSc linear solvers
    class PetscKSP
    {
    public:
        /// @brief Create a PetscKSP object
        PetscKSP();

        /// @brief Create a PetscMat from an existing PETSc KSP
        /// @param ksp Existing PETSc ksp
        /// @param inc_ref_count Whether to increase the ref count for ksp
        PetscKSP(KSP ksp, bool inc_ref_count);

        // Copy constructor (deleted)
        PetscKSP(const PetscKSP &) = delete;

        // Copy assignment operator (deleted)
        PetscKSP &operator=(PetscKSP &) = delete;

        /// @brief Move constructor
        PetscKSP(PetscKSP &&);

        /// @brief Move assignment
        PetscKSP &operator=(PetscKSP &&);

        /// @brief Destructor
        ~PetscKSP();

        /// @brief Get the underlying PETSc KSP
        KSP ksp() const;

        /// @brief Set the KSP options from the options database
        void set_from_options() const;

        /// @brief Set the prefix used for this specific KSP in the options database
        void set_options_prefix(const std::string &prefix) const;

        /// @brief Set the linear system LHS
        void set_operator(const PetscMat &A) const;

        /// @brief Solve the linear system Ax=b using the KSP
        int solve(const PetscVec &b, PetscVec &x) const;

    private:
        KSP ksp_;
    };
}

#endif // SFEM_HAS_PETSC