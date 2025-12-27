#ifdef SFEM_HAS_PETSC

#include <sfem/la/petsc/petsc_ksp.hpp>
#include <sfem/la/petsc/petsc_mat.hpp>
#include <sfem/la/petsc/petsc_vec.hpp>
#include <sfem/la/native/linear_solvers/linear_solver_factory.hpp>
#include <sfem/base/timer.hpp>
#include <sfem/base/error.hpp>

namespace sfem::la::petsc
{
    //=============================================================================
    PetscKSP::PetscKSP()
    {
        KSPCreate(PETSC_COMM_WORLD, &ksp_);
    }
    //=============================================================================
    PetscKSP::PetscKSP(KSP ksp, bool inc_ref_count)
        : ksp_(ksp)
    {
        if (!ksp_)
        {
            SFEM_ERROR("Invalid PETSc KSP passed to PetscKSP constructor");
        }

        if (inc_ref_count)
        {
            PetscObjectReference((PetscObject)ksp_);
        }
    }
    //=============================================================================
    PetscKSP::PetscKSP(PetscKSP &&other)
    {
        ksp_ = other.ksp_;
        other.ksp_ = nullptr;
    }
    //=============================================================================
    PetscKSP &PetscKSP::operator=(PetscKSP &&other)
    {
        if (this != &other)
        {
            if (ksp_)
            {
                KSPDestroy(&ksp_);
            }
            ksp_ = other.ksp_;
            other.ksp_ = nullptr;
        }

        return *this;
    }
    //=============================================================================
    PetscKSP::~PetscKSP()
    {
        if (ksp_)
        {
            KSPDestroy(&ksp_);
        }
    }
    //=============================================================================
    KSP PetscKSP::ksp() const
    {
        return ksp_;
    }
    //=============================================================================
    void PetscKSP::set_type(SolverType type) const
    {
        KSPType ksp_type;
        switch (type)
        {
        case SolverType::gmres:
            ksp_type = KSPGMRES;
            break;
        case SolverType::cg:
            ksp_type = KSPCG;
            break;
        default:
            ksp_type = KSPGMRES;
            break;
        }
        KSPSetType(ksp_, ksp_type);
    }
    //=============================================================================
    void PetscKSP::set_tolerances(SolverOptions options) const
    {
        KSPSetTolerances(ksp_,
                         options.rtol,
                         options.atol,
                         options.dtol,
                         options.n_iter_max);
    }
    //=============================================================================
    void PetscKSP::set_from_options() const
    {
        KSPSetFromOptions(ksp_);
    }
    //=============================================================================
    void PetscKSP::set_options_prefix(const std::string &prefix) const
    {
        KSPSetOptionsPrefix(ksp_, prefix.c_str());
    }
    //=============================================================================
    void PetscKSP::set_operator(const PetscMat &A) const
    {
        KSPSetOperators(ksp_, A.mat(), A.mat());
    }
    //=============================================================================
    bool PetscKSP::solve(const PetscVec &b, PetscVec &x) const
    {
        Timer timer;

        KSPType type;
        KSPGetType(ksp_, &type);
        if (strcmp(type, KSPPREONLY) != 0)
        {
            KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
        }

        KSPSetResidualHistory(ksp_, nullptr, PETSC_DECIDE, PETSC_TRUE);

        KSPSolve(ksp_, b.vec(), x.vec());

        // Get the number of iterations the KSP performed
        PetscInt n_iter;
        KSPGetIterationNumber(ksp_, &n_iter);

        // Check for convergence
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp_, &reason);
        if (reason < 0)
        {
            std::string msg = std::format("PetscKSP did not converge in {} iterations (reason {})\n",
                                          n_iter, static_cast<int>(reason));
            log_msg(msg, true, LogLevel::warning);
            return false;
        }
        else
        {
            return true;
        }
    }
    //=============================================================================
    std::vector<real_t> PetscKSP::residual_history() const
    {
        int n_iter;
        const real_t *hist;
        KSPGetResidualHistory(ksp_, &hist, &n_iter);

        std::vector<real_t> residual_history(n_iter);
        for (int i = 0; i < n_iter; i++)
        {
            residual_history[i] = hist[i];
        }

        return residual_history;
    }
}

#endif // SFEM_HAS_PETSC