#ifdef SFEM_HAS_PETSC

#include "petsc_linear_system.hpp"
#include "petsc.hpp"

namespace sfem::la::petsc
{
    //=============================================================================
    PetscLinearSystem::PetscLinearSystem(const IndexMap &index_map,
                                         const graph::Connectivity &connectivity,
                                         SolverType solver_type,
                                         SolverOptions solver_options,
                                         int block_size)
        : A_(create_mat(connectivity, index_map, index_map, block_size)),
          b_(create_vec(index_map, block_size)),
          solver_(),
          diag_(create_vec(index_map, block_size))
    {
        KSPType ksp_type;
        switch (solver_type)
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
        KSPSetType(solver_.ksp(), ksp_type);
        solver_.set_options(solver_options);
        solver_.set_from_options();
    }
    //=============================================================================
    PetscMat &PetscLinearSystem::A()
    {
        return A_;
    }
    //=============================================================================
    const PetscMat &PetscLinearSystem::A() const
    {
        return A_;
    }
    //=============================================================================
    PetscVec &PetscLinearSystem::b()
    {
        return b_;
    }
    //=============================================================================
    const PetscVec &PetscLinearSystem::b() const
    {
        return b_;
    }
    //=============================================================================
    PetscKSP &PetscLinearSystem::solver()
    {
        return solver_;
    }
    //=============================================================================
    const PetscKSP &PetscLinearSystem::solver() const
    {
        return solver_;
    }
    //=============================================================================
    void PetscLinearSystem::reset()
    {
        A_.zero_entries();
        b_.set_all(0.0);
    }
    //=============================================================================
    MatSet PetscLinearSystem::lhs()
    {
        return create_matset(A_);
    }
    //=============================================================================
    VecSet PetscLinearSystem::rhs()
    {
        return create_vecset(b_);
    }
    //=============================================================================
    void PetscLinearSystem::assemble()
    {
        A_.assemble();
        b_.assemble();
    }
    //=============================================================================
    void PetscLinearSystem::diagonal(la::Vector &diag) const
    {
        auto d = create_vec(diag);
        MatGetDiagonal(A_.mat(), d.vec());
    }
    //=============================================================================
    void PetscLinearSystem::scale_diagonal(real_t a)
    {
        MatGetDiagonal(A_.mat(), diag_.vec());
        VecScale(diag_.vec(), a);
        MatDiagonalSet(A_.mat(), diag_.vec(), INSERT_VALUES);
    }
    //=============================================================================
    void PetscLinearSystem::rhs_axpy(real_t a, const la::Vector &x)
    {
        auto x_ = create_vec(x);
        VecAXPY(b_.vec(), a, x_.vec());
    }
    //=============================================================================
    void PetscLinearSystem::eliminate_dofs(std::span<const int> idxs,
                                           std::span<const real_t> values)
    {
        diag_.set_values(idxs, values, INSERT_VALUES);
        eliminate_rows_cols(idxs, values, A_, b_, diag_);
    }
    //=============================================================================
    bool PetscLinearSystem::solve(Vector &x)
    {
        auto x_ = create_vec(x);
        solver_.set_operator(A_);
        return solver_.solve(b_, x_);
    }
}

#endif // SFEM_HAS_PETSC