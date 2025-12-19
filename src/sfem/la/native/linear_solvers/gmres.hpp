#pragma once

#include <sfem/la/native/linear_solvers/linear_solver.hpp>
#include <sfem/la/native/dense_matrix.hpp>
#include <sfem/la/native/vector.hpp>

namespace sfem::la
{
    /// @brief Generalized Minimum Residual solver
    class GMRES : public LinearSolver
    {
    public:
        GMRES(SolverOptions options = {}, int n_restart = 50);

    private:
        void init(const SparseMatrix &A, const Vector &b, Vector &x) override;

        void single_iteration(int iter, const SparseMatrix &A, const Vector &b, Vector &x) override;

        void restart(int iter, const SparseMatrix &A, const Vector &b, Vector &x);

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