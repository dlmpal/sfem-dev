#pragma once

#ifdef SFEM_HAS_SLEPC

#include "../petsc/petsc.hpp"
#include <slepc.h>

namespace sfem::la::slepc
{
    using namespace petsc;

    /// @brief Thin wrapper around the SLEPc's eigenproblem solvers (EPS)
    class SlepcEPS
    {
    public:
        /// @brief Create an EPS
        SlepcEPS();

        /// @brief Initialize from an existing EPS
        /// @param eps Existing EPS
        /// @param inc_ref_count Whether to increase the ref count the SLEPc object
        SlepcEPS(EPS eps, bool inc_ref_count);

        // Avoid copies
        SlepcEPS(const SlepcEPS &) = delete;
        SlepcEPS &operator=(SlepcEPS &) = delete;

        // Move constructor and assignment
        SlepcEPS(SlepcEPS &&);
        SlepcEPS &operator=(SlepcEPS &&);

        // Destructor
        ~SlepcEPS();

        /// @brief Get the underlying SLEPc EPS
        EPS eps() const;

        /// @brief Set the EPS options from the options database
        void set_from_options() const;

        /// @brief Set the operators for the generalized eigenproblem
        void set_operators(const PetscMat &A,
                           const PetscMat &B) const;

        /// @brief Set the operator for the standard eigenproblem
        void set_operators(const PetscMat &A) const;

        /// @brief Solve the eigenvalue problem using the EPS
        /// @param n_pairs Number of eigenpairs to compute
        int solve(int n_pairs) const;

        /// @brief Get the number of converged eigenpairs
        int n_converged() const;

        /// @brief Get the i-th eigenvalue (real and imaginary part)
        std::array<real_t, 2> eigenvalue(int pair_idx) const;

        /// @brief Get the i-th eigenpair
        /// @return Real and imaginary parts of eigenvalue and eigenvector
        std::pair<std::array<real_t, 2>, std::array<PetscVec, 2>>
        eigenpair(int pair_idx) const;

    private:
        EPS eps_;
    };
}

#endif // SFEM_HAS_SLEPC