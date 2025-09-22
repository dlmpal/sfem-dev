#pragma once

#include "linear_solver.hpp"
#include "../vector.hpp"
#include "../dense_matrix.hpp"

namespace sfem::la
{

    class GMRES : public LinearSolver
    {
    public:
        GMRES(real_t tol, int n_iter_max, bool verbose, int n_restart);

    private:
        void init(const SparseMatrix &A,
                  const Vector &b, Vector &x) override;

        void single_iteration(int iter, const SparseMatrix &A,
                              const Vector &b, Vector &x) override;

        void init_(int iter, const SparseMatrix &A,
                   const Vector &b, Vector &x);

    private:
        /// @brief Number of iterations before restart
        int n_restart_;

        /// @brief Iteration since last restart
        int riter_;

        // Workspace vectors
        Vector v1;
        Vector v2;
        Vector x0;

        DenseMatrix Q;
        DenseMatrix H;
        DenseMatrix e1;
    };
}