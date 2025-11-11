#pragma once

#include "linear_solver.hpp"
#include "../vector.hpp"

namespace sfem::la
{
    /// @brief Conjugate Gradient solver
    class CG : public LinearSolver
    {
    public:
        CG(const SolverOptions &options);

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