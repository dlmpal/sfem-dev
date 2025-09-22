#pragma once

#include "linear_solver.hpp"
#include "../vector.hpp"

namespace sfem::la
{
    class CG : public LinearSolver
    {
    public:
        CG(real_t tol, int n_iter_max, bool verbose);

    private:
        void single_iteration(int iter, const SparseMatrix &A,
                              const Vector &b, Vector &x) override;

        void init(const SparseMatrix &A,
                  const Vector &b, Vector &x) override;

    private:
        /// @brief  Workspace vector, used for storing intermediate products
        Vector Ap;

        // Search direction vector
        Vector p;

        // Residual vector
        Vector r;
    };
}