#pragma once

#include "../../../base/config.hpp"
#include <string>
#include <vector>

namespace sfem::la
{
    // Forward declarations
    class Vector;
    class SparseMatrix;

    /// @brief Linear solver ABC
    class LinearSolver
    {
    public:
        /// @brief Create a solver
        LinearSolver(const std::string &name, real_t tol,
                     int n_iter_max, bool verbose);

        /// @brief Get the solver's name
        std::string name() const;

        /// @brief Get the solver's termination tolernace
        real_t tol() const;

        /// @brief Get the solver's maximum allowed number of iterations
        int n_iter_max() const;

        /// @brief Get the solver's verbosity
        bool verbose() const;

        /// @brief Get the solver's residual history
        std::vector<real_t> residual_history() const;

        /// @brief Run the solver, i.e. solve Ax=b for x
        void run(const SparseMatrix &A,
                 const Vector &b, Vector &x);

    protected:
        /// @brief Initialize various solver attributes such as workspace vectors.
        /// @note Should also compute the first (0-th) residual
        virtual void init(const SparseMatrix &A,
                          const Vector &b, Vector &x) = 0;

        /// @brief Perform a single solver iteration, updating the solution values
        /// @note Should also update residual history
        virtual void single_iteration(int iter, const SparseMatrix &A,
                                      const Vector &b, Vector &x) = 0;

    protected:
        /// @brief Linear solver name
        std::string name_;

        /// @brief Termination tolerance
        real_t tol_;

        /// @brief Number of maximum iterations
        int n_iter_max_;

        /// @brief Verbosity
        bool verbose_;

        /// @brief Residual norm history
        mutable std::vector<real_t> residual_history_;
    };
}