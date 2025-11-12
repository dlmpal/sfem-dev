#include "linear_solver.hpp"
#include "../../../base/error.hpp"
#include "../vector.hpp"
#include "../sparse_matrix.hpp"

namespace sfem::la
{
    //=============================================================================
    LinearSolver::LinearSolver(const std::string &name, SolverOptions options)
        : name_(name),
          options_(options)
    {
    }
    //=============================================================================
    std::string LinearSolver::name() const
    {
        return name_;
    }
    //=============================================================================
    SolverOptions &LinearSolver::options()
    {
        return options_;
    }
    //======================================================================
    SolverOptions LinearSolver::options() const
    {
        return options_;
    }
    //=============================================================================
    std::vector<real_t> LinearSolver::residual_history() const
    {
        return residual_history_;
    }
    //=============================================================================
    bool LinearSolver::run(const SparseMatrix &A, const Vector &b, Vector &x)
    {
        // Check that options are valid
        if (options_.atol < 0)
        {
            SFEM_ERROR(std::format("Invalid absolute tolerance {} (<0)\n", options_.atol));
        }
        if (options_.rtol < 0)
        {
            SFEM_ERROR(std::format("Invalid relative tolerance {} (<0)\n", options_.rtol));
        }
        if (options_.dtol < 0)
        {
            SFEM_ERROR(std::format("Invalid divergence tolerance {} (<0)\n", options_.dtol));
        }
        if (options_.n_iter_max <= 0)
        {
            SFEM_ERROR(std::format("Invalid number of iterations {} (<=0)\n", options_.n_iter_max));
        }

        // Reset residual history
        residual_history_.resize(options_.n_iter_max + 1, 0.0);

        // Initialize the solver
        init(A, b, x);
        if (options_.print_iter)
        {
            log_msg(std::format("{} - Iteration 0, Residual {}\n",
                                name_, residual_history_[0]),
                    true);
        }

        // Compute the termination tolerance
        const real_t r0 = residual_history_[0];
        const real_t tol = std::max(options_.atol, options_.rtol * r0);

        // Perform iterations
        int iter = 0;
        while (residual_history_[iter] >= tol && iter < options_.n_iter_max)
        {
            iter++;
            single_iteration(iter, A, b, x);

            if (options_.print_iter)
            {
                log_msg(std::format("{} Iteration {}, Residual {}\n",
                                    name_, iter, residual_history_[iter]),
                        true);
            }

            // Check for divergence
            if (residual_history_[iter] >= options_.dtol * r0)
            {
                log_msg(std::format("{} has diverged in {} iterations\n", name_, iter), true);
                return false;
            }
        }

        // Check for convergence
        bool converged = residual_history_[iter] < tol ? true : false;

        // Print convergence message
        if (options_.print_conv)
        {
            log_msg(std::format("{} Initial Residual {}, Final Residual {}\n", name_, r0, residual_history_[iter]), true);

            if (converged)
            {
                log_msg(std::format("{} has converged in {} iterations\n", name_, iter), true);
            }
            else
            {
                log_msg(std::format("{} has failed to converge in {} iterations. Residual ({}) is greater than tolerance ({})\n",
                                    name_, options_.n_iter_max, residual_history_[iter], tol),
                        true);
            }
        }

        return converged;
    }
}