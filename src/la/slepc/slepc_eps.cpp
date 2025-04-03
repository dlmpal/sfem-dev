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
        SFEM_CHECK_PETSC_ERROR(EPSCreate(MPI_COMM_WORLD, &eps_));
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
            SFEM_CHECK_PETSC_ERROR(PetscObjectReference((PetscObject)eps_));
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
                SFEM_CHECK_PETSC_ERROR(EPSDestroy(&eps_));
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
            SFEM_CHECK_PETSC_ERROR(EPSDestroy(&eps_));
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
        SFEM_CHECK_PETSC_ERROR(EPSSetFromOptions(eps_));
    }
    //=============================================================================
    void SlepcEPS::set_operators(const PetscMat &A, const PetscMat &B) const
    {
        SFEM_CHECK_PETSC_ERROR(EPSSetOperators(eps_, A.mat(), B.mat()));
        SFEM_CHECK_PETSC_ERROR(EPSSetProblemType(eps_, EPS_GHEP));
    }
    //=============================================================================
    void SlepcEPS::set_operators(const PetscMat &A) const
    {
        SFEM_CHECK_PETSC_ERROR(EPSSetOperators(eps_, A.mat(), nullptr));
        SFEM_CHECK_PETSC_ERROR(EPSSetProblemType(eps_, EPS_HEP));
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
        SFEM_CHECK_PETSC_ERROR(EPSSetDimensions(eps_, n_pairs,
                                                PETSC_DECIDE, PETSC_DECIDE));

        // Solve eigenproblem
        SFEM_CHECK_PETSC_ERROR(EPSSolve(eps_));

        // Get the number of iterations the EPS performed
        int n_iter;
        SFEM_CHECK_PETSC_ERROR(EPSGetIterationNumber(eps_, &n_iter));

        // Check for convergence
        EPSConvergedReason reason;
        SFEM_CHECK_PETSC_ERROR(EPSGetConvergedReason(eps_, &reason));
        if (reason < 0)
        {
            std::string msg = std::format("SLEPcEPS did not converge in {} iterations (reason {})\n",
                                          n_iter, static_cast<int>(reason));
            log_msg(msg, LogLevel::warning);
        }

        return static_cast<int>(n_iter);
    }
    //=============================================================================
    int SlepcEPS::n_converged() const
    {
        int n_pairs;
        SFEM_CHECK_PETSC_ERROR(EPSGetConverged(eps_, &n_pairs));
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
        int error_code = EPSGetEigenpair(eps_, pair_idx,
                                         &real, &imag,
                                         nullptr, nullptr);
        SFEM_CHECK_PETSC_ERROR(error_code);

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
        SFEM_CHECK_PETSC_ERROR(EPSGetOperators(eps_, &mat, nullptr));
        SFEM_CHECK_PETSC_ERROR(MatGetSize(mat, &n_global, &n_local));

        // Get the eigenvalue and its corresponding eigenvector
        real_t real, imag;
        auto V_real = petsc::create_vec(IndexMap(n_local), 1);
        auto V_imag = V_real.copy();
        int error_code = EPSGetEigenpair(eps_, pair_idx,
                                         &real, &imag,
                                         V_real.vec(), V_imag.vec());
        SFEM_CHECK_PETSC_ERROR(error_code);

        return {std::array{real, imag},
                std::array{std::move(V_real), std::move(V_imag)}};
    }
}

#endif // SFEM_HAS_SLEPC