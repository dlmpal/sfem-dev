#ifdef SFEM_HAS_SLEPC

#include "slepc_eps.hpp"
#include "../petsc/petsc_vec.hpp"
#include "../petsc/petsc_mat.hpp"
#include "../../base/error.hpp"
#include "../../base/timer.hpp"

namespace sfem::la::slepc
{
    //=============================================================================
    SlepcEPS::SlepcEPS()
    {
        EPSCreate(MPI_COMM_WORLD, &eps_);
    }
    //=============================================================================
    SlepcEPS::SlepcEPS(EPS eps, bool inc_ref_count)
        : eps_(eps)
    {
        if (!eps_)
        {
            SFEM_ERROR("Invalid SLEPc EPS passed to SlepcEPS constructor");
        }

        if (inc_ref_count)
        {
            PetscObjectReference((PetscObject)eps_);
        }
    }
    //=============================================================================
    SlepcEPS::SlepcEPS(SlepcEPS &&other)
    {
        eps_ = other.eps_;
        other.eps_ = nullptr;
    }
    //=============================================================================
    SlepcEPS &SlepcEPS::operator=(SlepcEPS &&other)
    {
        if (this != &other)
        {
            if (eps_)
            {
                EPSDestroy(&eps_);
            }
            eps_ = other.eps_;
            other.eps_ = nullptr;
        }

        return *this;
    }
    //=============================================================================
    SlepcEPS::~SlepcEPS()
    {
        if (eps_)
        {
            EPSDestroy(&eps_);
        }
    }
    //=============================================================================
    EPS SlepcEPS::eps() const
    {
        return eps_;
    }
    //=============================================================================
    void SlepcEPS::set_from_options() const
    {
        EPSSetFromOptions(eps_);
    }
    //=============================================================================
    void SlepcEPS::set_operators(const PetscMat &A, const PetscMat &B) const
    {
        EPSSetOperators(eps_, A.mat(), B.mat());
        EPSSetProblemType(eps_, EPS_GHEP);
    }
    //=============================================================================
    void SlepcEPS::set_operators(const PetscMat &A) const
    {
        EPSSetOperators(eps_, A.mat(), nullptr);
        EPSSetProblemType(eps_, EPS_HEP);
    }
    //=============================================================================
    int SlepcEPS::solve(int n_pairs) const
    {
        Timer timer;

        // Set EPS dimensions
        if (n_pairs <= 0)
        {
            SFEM_ERROR(std::format("Invalid number of pairs {} (<=0)\n", n_pairs));
        }
        EPSSetDimensions(eps_, n_pairs, PETSC_DECIDE, PETSC_DECIDE);

        // Solve eigenproblem
        EPSSolve(eps_);

        // Get the number of iterations the EPS performed
        int n_iter;
        EPSGetIterationNumber(eps_, &n_iter);

        // Check for convergence
        EPSConvergedReason reason;
        EPSGetConvergedReason(eps_, &reason);
        if (reason < 0)
        {
            std::string msg = std::format("SLEPcEPS did not converge in {} iterations (reason {})\n",
                                          n_iter, static_cast<int>(reason));
            log_msg(msg, true, LogLevel::warning);
        }

        return static_cast<int>(n_iter);
    }
    //=============================================================================
    int SlepcEPS::n_converged() const
    {
        int n_pairs;
        EPSGetConverged(eps_, &n_pairs);
        return static_cast<int>(n_pairs);
    }
    //=============================================================================
    std::array<real_t, 2> SlepcEPS::eigenvalue(int pair_idx) const
    {
        // Check if the pair index is greater than the number of converged pairs
        if (pair_idx >= n_converged())
        {
            SFEM_ERROR(std::format("Eigenvalue {} is not available\n", pair_idx));
        }

        // Get the eigenvalue
        real_t real, imag;
        EPSGetEigenpair(eps_, pair_idx, &real, &imag, nullptr, nullptr);

        return {real, imag};
    }
    //=============================================================================
    std::pair<std::array<real_t, 2>, std::array<petsc::PetscVec, 2>>
    SlepcEPS::eigenpair(int pair_idx) const
    {
        // Check if the pair index is greater than the number of converged pairs
        if (pair_idx >= n_converged())
        {
            SFEM_ERROR(std::format("Eigenpair {} is not available\n", pair_idx));
        }

        // Get operator local and global sizes
        Mat mat;
        int n_local, n_global;
        EPSGetOperators(eps_, &mat, nullptr);
        MatGetSize(mat, &n_global, &n_local);

        // Get the eigenvalue and its corresponding eigenvector
        real_t real, imag;
        auto V_real = petsc::create_vec(IndexMap(n_local), 1);
        auto V_imag = V_real.copy();
        EPSGetEigenpair(eps_, pair_idx,
                        &real, &imag,
                        V_real.vec(), V_imag.vec());

        return {std::array{real, imag},
                std::array{std::move(V_real), std::move(V_imag)}};
    }
}

#endif // SFEM_HAS_SLEPC