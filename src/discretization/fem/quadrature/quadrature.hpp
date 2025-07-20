#pragma once

#include "../../../base/config.hpp"
#include <array>

namespace sfem::fem
{
    // Forward declaration
    class IntegrationPoint;

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

        /// @brief Get the integration weight and point coordinates
        virtual IntegrationPoint point(int i) const = 0;

    protected:
        /// @brief Number of integration points
        int n_points_;
    };

    /// @brief Integration point defined by
    /// point's weight and coordinates
    struct IntegrationPoint
    {
        real_t x() const
        {
            return point[0];
        }

        real_t y() const
        {
            return point[1];
        }

        real_t z() const
        {
            return point[2];
        }

        /// @brief Integration point weight
        real_t weight{};

        /// @brief Integration point coordinates
        std::array<real_t, 3> point{};
    };
}