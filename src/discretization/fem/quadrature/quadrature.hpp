#pragma once

#include "../../../base/config.hpp"
#include <array>

namespace sfem::fem
{
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
        std::array<real_t, 3> point;
    };

    /// @brief Integration rule ABC
    class IntegrationRule
    {
    public:
        IntegrationRule(int order)
            : order_(order)
        {
        }

        int order() const
        {
            return order_;
        }

        int &order()
        {
            return order_;
        }

        virtual ~IntegrationRule() = default;
        virtual int n_points() const = 0;
        virtual IntegrationPoint point(int i) const = 0;

    protected:
        int order_;
    };
}