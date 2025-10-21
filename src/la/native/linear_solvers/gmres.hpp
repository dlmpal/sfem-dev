#pragma once

#include "linear_solver.hpp"
#include "../vector.hpp"
#include "../dense_matrix.hpp"

namespace sfem::la
{

    class GMRES : public LinearSolver
    {
    public:
        GMRES(real_t tol, int n_iter_max, bool verbose, int n_restart = 50);

    private:
        void init(const SparseMatrix &A,
                  const Vector &b, Vector &x) override;

        void single_iteration(int iter, const SparseMatrix &A,
                              const Vector &b, Vector &x) override;

        void restart(int iter, const SparseMatrix &A,
                     const Vector &b, Vector &x);

    private:
        /// @brief Number of iterations before restart
        int n_restart_;

        /// @brief Iterations since last restart
        int riter_;

        /// @brief Initial solution vector (since last restart)
        Vector x0_;

        /// @brief Krylov subspace orthonormal basis vectors
        std::vector<Vector> Q_;

        /// @brief Hessenberg matrix
        DenseMatrix H_;

        /// @brief e1 vector (e1=[||r0||, 0, 0, ..., 0]^T)
        DenseMatrix e1_;
    };
}