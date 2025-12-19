#pragma once

#include <sfem/la/native/linear_solvers/linear_solver.hpp>
#include <sfem/la/native/vector.hpp>

namespace sfem::la
{
    /// @brief Conjugate Gradient solver
    class CG : public LinearSolver
    {
    public:
        CG(SolverOptions options = {});

    private:
        void init(const SparseMatrix &A, const Vector &b, Vector &x) override;

        void single_iteration(int iter, const SparseMatrix &A, const Vector &b, Vector &x) override;

    private:
        /// @brief  Workspace vector, used for storing intermediate products
        Vector Ap;

        // Search direction vector
        Vector p;

        // Residual vector
        Vector r;
    };
}