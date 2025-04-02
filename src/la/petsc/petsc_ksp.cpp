#ifdef SFEM_HAS_PETSC

#include "petsc_ksp.hpp"
#include "petsc_mat.hpp"
#include "petsc_vec.hpp"
#include "../../base/timer.hpp"
#include "../../base/error.hpp"

namespace sfem::la::petsc
{
    //=============================================================================
    PetscKSP::PetscKSP()
    {
        SFEM_CHECK_PETSC_ERROR(KSPCreate(PETSC_COMM_WORLD, &ksp_));
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
            SFEM_CHECK_PETSC_ERROR(PetscObjectReference((PetscObject)ksp_));
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
                SFEM_CHECK_PETSC_ERROR(KSPDestroy(&ksp_));
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
            SFEM_CHECK_PETSC_ERROR(KSPDestroy(&ksp_));
        }
    }
    //=============================================================================
    KSP PetscKSP::ksp() const
    {
        return ksp_;
    }
    //=============================================================================
    void PetscKSP::set_from_options() const
    {
        SFEM_CHECK_PETSC_ERROR(KSPSetFromOptions(ksp_));
    }
    //=============================================================================
    void PetscKSP::set_options_prefix(const std::string &prefix) const
    {
        SFEM_CHECK_PETSC_ERROR(KSPSetOptionsPrefix(ksp_, prefix.c_str()));
    }
    //=============================================================================
    void PetscKSP::set_operator(const PetscMat &A) const
    {
        SFEM_CHECK_PETSC_ERROR(KSPSetOperators(ksp_, A.mat(), A.mat()));
    }
    //=============================================================================
    int PetscKSP::solve(const PetscVec &b, PetscVec &x) const
    {
        Timer timer;

        SFEM_CHECK_PETSC_ERROR(KSPSolve(ksp_, b.vec(), x.vec()));

        // Get the number of iterations the KSP performed
        PetscInt n_iter;
        SFEM_CHECK_PETSC_ERROR(KSPGetIterationNumber(ksp_, &n_iter));

        // Check for convergence
        KSPConvergedReason reason;
        SFEM_CHECK_PETSC_ERROR(KSPGetConvergedReason(ksp_, &reason));
        if (reason < 0)
        {
            std::string msg = std::format("PetscKSP did not converge in {} iterations (reason {})\n",
                                          n_iter, static_cast<int>(reason));
            log_msg(msg, LogLevel::warning);
        }

        return static_cast<int>(n_iter);
    }
}

#endif // SFEM_HAS_PETSC