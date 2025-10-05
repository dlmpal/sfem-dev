#include "erk.hpp"
#include "../../la/native/vector.hpp"
#include "../../base/error.hpp"

namespace sfem::ode
{
    //=============================================================================
    ERKIntegrator::ERKIntegrator(const la::Vector &state, RHSFunction rhs,
                                 int n_stages, std::vector<real_t> &&nodes,
                                 std::vector<real_t> &&weights,
                                 la::DenseMatrix &&coeffs)
        : rhs_(rhs),
          n_stages_(n_stages),
          nodes_(std::move(nodes)),
          weights_(std::move(weights)),
          coeffs_(std::move(coeffs))
    {
        SFEM_CHECK_SIZES(nodes_.size(), weights_.size());
        SFEM_CHECK_SIZES(nodes_.size(), coeffs.n_rows());
        SFEM_CHECK_SIZES(nodes_.size(), coeffs.n_cols());

        for (int i = 0; i < n_stages_; i++)
        {
            stages_.emplace_back(state.index_map(), state.block_size());
        }
    }
    //=============================================================================
    void ERKIntegrator::advance(const la::Vector &S_old, la::Vector &S_new, real_t time, real_t dt)
    {
        // Evaluate the RHS for each stage
        for (int i = 0; i < n_stages_; i++)
        {
            // Evaluate the intermediate state
            la::copy(S_old, S_new);
            for (int j = 0; j < i; j++)
            {
                real_t aij = coeffs_(i, j);
                la::axpy(dt * aij, stages_[j], S_new);
            }
            S_new.update_ghosts();
            rhs_(S_new, stages_[i], time + dt * nodes_[i]);
        }

        // Evaluate the next state
        la::copy(S_old, S_new);
        for (int i = 0; i < n_stages_; i++)
        {
            la::axpy(dt * weights_[i], stages_[i], S_new);
        }
    }
    //=============================================================================
    static ERKIntegrator create_fe(const la::Vector &state, ERKIntegrator::RHSFunction rhs)
    {
        const int n_stages = 1;
        std::vector<real_t> nodes = {0.};
        std::vector<real_t> weights = {1.0};
        la::DenseMatrix coeffs(n_stages, n_stages, 0.);

        return ERKIntegrator(state, rhs, n_stages,
                             std::move(nodes),
                             std::move(weights),
                             std::move(coeffs));
    }
    //=============================================================================
    static ERKIntegrator create_rk4(const la::Vector &state, ERKIntegrator::RHSFunction rhs)
    {
        const int n_stages = 4;
        std::vector<real_t> nodes = {0., 1. / 2., 1. / 2., 1.};
        std::vector<real_t> weights = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};
        la::DenseMatrix coeffs(n_stages, n_stages, 0);
        coeffs(1, 0) = 0.5;
        coeffs(2, 1) = 0.5;
        coeffs(3, 2) = 1.0;

        return ERKIntegrator(state, rhs, n_stages,
                             std::move(nodes),
                             std::move(weights),
                             std::move(coeffs));
    }
    //=============================================================================
    ERKIntegrator create_erk(const la::Vector &state, ERKIntegrator::RHSFunction rhs, ERKType type)
    {
        switch (type)
        {
        case ERKType::fe:
            return create_fe(state, rhs);
        case ERKType::rk4:
            return create_rk4(state, rhs);
        default:
            SFEM_ERROR("Invalid ERK type\n");
            return ERKIntegrator(state, rhs, 0, {}, {}, la::DenseMatrix(1, 1));
        }
    }
}