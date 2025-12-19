#pragma once

#include <sfem/base/config.hpp>
#include <string>
#include <vector>

namespace sfem::la
{
    // Forward declarations
    class Vector;
    class SparseMatrix;

    /// @brief Linear solver options
    struct SolverOptions
    {
        /// @brief Absolute tolerance
        real_t atol = 1e-10;

        /// @brief Relative tolerance
        real_t rtol = 1e-6;

        /// @brief Divergence tolerance
        real_t dtol = 1e3;

        /// @brief Number of maximum iterations
        int n_iter_max = 500;

        /// @brief Whether to print convergence related info
        bool print_conv = true;

        /// @brief Whether to print iteration info
        bool print_iter = false;
    };

    /// @brief Linear solver ABC
    class LinearSolver
    {
    public:
        /// @brief Create a solver
        LinearSolver(const std::string &name, SolverOptions options);

        /// @brief Get the solver's name
        std::string name() const;

        /// @brief Get the solver's options
        SolverOptions &options();

        /// @brief Get the solver's options (const version)
        SolverOptions options() const;

        /// @brief Get the solver's residual history
        std::vector<real_t> residual_history() const;

        /// @brief Run the solver, i.e. solve Ax=b for x
        bool run(const SparseMatrix &A, const Vector &b, Vector &x);

    protected:
        /// @brief Initialize various solver attributes such as workspace vectors.
        /// @note Should also compute the first (0-th) residual
        virtual void init(const SparseMatrix &A, const Vector &b, Vector &x) = 0;

        /// @brief Perform a single solver iteration, updating the solution values
        /// @note Should also update residual history
        virtual void single_iteration(int iter, const SparseMatrix &A, const Vector &b, Vector &x) = 0;

    protected:
        /// @brief Solver name
        std::string name_;

        /// @brief Solver options
        SolverOptions options_;

        /// @brief Residual norm history
        std::vector<real_t> residual_history_;
    };
}