#include "linear_solver.hpp"
#include "../../../base/error.hpp"

namespace sfem::la
{
    //=============================================================================
    LinearSolver::LinearSolver(const std::string &name, real_t tol, int n_iter_max, bool verbose)
        : name_(name),
          tol_(tol),
          n_iter_max_(n_iter_max),
          verbose_(verbose),
          residual_history_(n_iter_max_)
    {
        /// @todo check values of tole and n_iter_max
    }
    //=============================================================================
    std::string LinearSolver::name() const
    {
        return name_;
    }
    //=============================================================================
    real_t LinearSolver::tol() const
    {
        return tol_;
    }
    //=============================================================================

    int LinearSolver::n_iter_max() const
    {
        return n_iter_max_;
    }
    //=============================================================================
    bool LinearSolver::verbose() const
    {
        return verbose_;
    }
    //=============================================================================
    std::vector<real_t> LinearSolver::residual_history() const
    {
        return residual_history_;
    }
    //=============================================================================
    void LinearSolver::run(const SparseMatrix &A,
                           const Vector &b,
                           Vector &x)
    {
        // Reset residual history
        std::fill(residual_history_.begin(),
                  residual_history_.end(), 0.0);

        // Initialize the solver
        init(A, b, x);
        if (verbose_)
        {
            log_msg(std::format("{} - Iteration 0, Residual {}\n",
                                name_, residual_history_[0]),
                    true);
        }

        int iter = 1;
        bool converged = false;
        for (; iter < n_iter_max_; iter++)
        {
            // Perform a single iteration
            single_iteration(iter, A, b, x);

            // Print current iteration number and residual
            if (verbose_)
            {
                log_msg(std::format("{} - Iteration {}, Residual {}\n",
                                    name_, iter, residual_history_[iter]),
                        true);
            }

            // Check for convergence
            if (residual_history_[iter] <= tol_)
            {
                converged = true;
                break;
            }
        }

        // Print convergence message
        if (verbose_)
        {
            if (converged)
            {
                log_msg(std::format("{} has converged in {} iterations\n",
                                    name_, iter - 1),
                        true);
            }
            else
            {
                log_msg(std::format("{} has failed to converge in {} iterations. Residual ({}) is greater than tolerance ({})\n",
                                    name_, n_iter_max_, residual_history_[iter - 1], tol_),
                        true);
            }
        }
    }
}