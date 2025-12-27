#include "linear_system.hpp"

namespace sfem::la
{
    //=============================================================================
    NativeLinearSystem::NativeLinearSystem(std::shared_ptr<const IndexMap> index_map,
                                           std::shared_ptr<const graph::Connectivity> connectivity,
                                           SolverType solver_type,
                                           SolverOptions solver_options,
                                           int block_size)
        : A_(connectivity, index_map, index_map, block_size),
          b_(index_map, block_size),
          solver_(create_solver(solver_type, solver_options))
    {
    }
    //=============================================================================
    SparseMatrix &NativeLinearSystem::A()
    {
        return A_;
    }
    //=============================================================================
    const SparseMatrix &NativeLinearSystem::A() const
    {
        return A_;
    }
    //=============================================================================
    Vector &NativeLinearSystem::b()
    {
        return b_;
    }
    //=============================================================================
    const Vector &NativeLinearSystem::b() const
    {
        return b_;
    }
    //=============================================================================
    std::shared_ptr<LinearSolver> NativeLinearSystem::solver()
    {
        return solver_;
    }
    //=============================================================================
    std::shared_ptr<const LinearSolver> NativeLinearSystem::solver() const
    {
        return solver_;
    }
    //=============================================================================
    void NativeLinearSystem::reset()
    {
        A_.set_all(0.0);
        b_.set_all(0.0);
    }
    //=============================================================================
    MatSet NativeLinearSystem::lhs()
    {
        return create_matset(A_);
    }
    //=============================================================================
    VecSet NativeLinearSystem::rhs()
    {
        return create_vecset(b_);
    }
    //=============================================================================
    void NativeLinearSystem::assemble()
    {
        A_.assemble();
        b_.assemble();
    }
    //=============================================================================
    void NativeLinearSystem::diagonal(la::Vector &diag) const
    {
        A_.diagonal(diag);
    }
    //=============================================================================
    void NativeLinearSystem::scale_diagonal(real_t a)
    {
        A_.scale_diagonal(a);
    }
    //=============================================================================
    void NativeLinearSystem::rhs_axpy(real_t a, const la::Vector &x)
    {
        axpy(a, x, b_);
    }
    //=============================================================================
    void NativeLinearSystem::eliminate_dofs(std::span<const int>,
                                            std::span<const real_t>)
    {
        /// @todo
    }
    //=============================================================================
    bool NativeLinearSystem::solve(Vector &x)
    {
        return solver_->run(A_, b_, x);
    }
    //=============================================================================
    std::vector<real_t> NativeLinearSystem::residual_history() const
    {
        return solver_->residual_history();
    }
}