#pragma once

#include "../../../base/config.hpp"
#include <array>

namespace sfem::fem
{
    // Forward declaration
    // class IntegrationPoint;

    /// @brief Integration rule ABC
    class IntegrationRule
    {
    public:
        /// @brief Create a rule for a given number of points
        IntegrationRule(int n_points);

        // Destructor
        virtual ~IntegrationRule() = default;

        /// @brief Get the number of integration points
        int n_points() const;

        /// @brief Set the number of integration points
        /// @note For Gauss quadrature this sets the no. points
        /// per direction
        virtual void set_n_points(int n_points);

        /// @brief Get the integration point weight
        virtual real_t weight(int qpt_idx) const = 0;

        /// @brief Get the integration point coordinates
        virtual std::array<real_t, 3> point(int qpt_idx) const = 0;

    protected:
        /// @brief Number of integration points
        int n_points_;
    };
}