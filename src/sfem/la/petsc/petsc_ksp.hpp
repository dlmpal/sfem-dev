#pragma once

#ifdef SFEM_HAS_PETSC

#include <sfem/la/petsc/petsc_vec.hpp>
#include <sfem/la/petsc/petsc_mat.hpp>
#include <petscksp.h>

// Forward declaration
namespace sfem::la
{
    enum class SolverType;
    struct SolverOptions;
}

namespace sfem::la::petsc
{
    /// @brief Thin wrapper around PETSc's linear solvers (KSP)
    class PetscKSP
    {
    public:
        /// @brief Create a KSP
        PetscKSP();

        /// @brief Initialize from an existing KSP
        /// @param ksp Existing KSP
        /// @param inc_ref_count Whether to increase the ref count for the PETSc object
        PetscKSP(KSP ksp, bool inc_ref_count);

        // Avoid copies
        PetscKSP(const PetscKSP &) = delete;
        PetscKSP &operator=(PetscKSP &) = delete;

        // Move constructor and assignment
        PetscKSP(PetscKSP &&);
        PetscKSP &operator=(PetscKSP &&);

        // Destructor
        ~PetscKSP();

        /// @brief Get the underlying KSP
        KSP ksp() const;

        /// @brief Set the KSP tolerances using native solver types
        /// @note See SolverType for currently available native types
        void set_type(SolverType type) const;

        /// @brief Set the KSP tolerances using native solver options
        /// @note Not all entries in options are utilized
        void set_tolerances(SolverOptions options) const;

        /// @brief Set the KSP options from the options database
        void set_from_options() const;

        /// @brief Set the prefix used for this specific KSP in the options database
        void set_options_prefix(const std::string &prefix) const;

        /// @brief Set the linear system LHS
        void set_operator(const PetscMat &A) const;

        /// @brief Solve the linear system Ax=b using the KSP
        bool solve(const PetscVec &b, PetscVec &x) const;

        std::vector<real_t> residual_history() const;

    private:
        /// @brief Underlying PETSc KSP
        KSP ksp_;
    };
}

#endif // SFEM_HAS_PETSC