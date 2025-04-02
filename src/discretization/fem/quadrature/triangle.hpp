#pragma once

#include "quadrature.hpp"

namespace sfem::fem::quadrature
{
    class Triangle : public IntegrationRule
    {
    public:
        Triangle(int order)
            : IntegrationRule(order)
        {
        }
        int n_points() const override;
        IntegrationPoint point(int i) const override;
    };
}