#pragma once

#include "../../la/native/dense_matrix.hpp"
#include <functional>

// Forward declaration
namespace sfem::la
{
    class Vector;
}

namespace sfem::ode
{
    /// @brief Explicit Runge-Kutta integrator
    class ERKIntegrator
    {
    public:
        /// @brief Right-hand-side function F = F(t, S)
        using RHSFunction = std::function<void(const la::Vector &S, la::Vector &rhs, real_t time)>;

        /// @brief Create an Explicit Runge-Kutta integrator
        /// @param state State vector
        /// @param rhs RHS function
        /// @param n_stages Number of stages
        /// @param nodes Butcher tableau nodes
        /// @param weights Butcher tableau weights
        /// @param coeffs Butcher tableau coefficients
        ERKIntegrator(const la::Vector &state, RHSFunction rhs,
                      int n_stages, std::vector<real_t> &&nodes,
                      std::vector<real_t> &&weights,
                      la::DenseMatrix &&coeffs);

        /// @brief Advance the state forward by one timestep
        /// @param S_old Old timestep state vector
        /// @param S_new New timestep state vector
        /// @param time Current time
        /// @param dt Current timestep size
        void advance(const la::Vector &S_old, la::Vector &S_new, real_t time, real_t dt);

    private:
        /// @brief RHS function
        RHSFunction rhs_;

        /// @brief Number of stages
        int n_stages_;

        /// @brief Butcher tableau nodes (ci)
        std::vector<real_t> nodes_;

        /// @brief Butcher tableau weights (bi)
        std::vector<real_t> weights_;

        /// @brief Butcher tableau coefficients (aij)
        la::DenseMatrix coeffs_;

        /// @brief RHS vector for the intermediate stages
        std::vector<la::Vector> stages_;
    };

    /// @brief Available Explicit Runge-Kutta integrator types
    enum class ERKType
    {
        fe, // Forward Euler
        rk4 // Classic fourth-order Runge-Kutta
    };

    /// @brief Create an ERK integrator of specific type
    ERKIntegrator create_erk(const la::Vector &state, ERKIntegrator::RHSFunction rhs, ERKType type);
}