#pragma once

#include "quadrature.hpp"

namespace sfem::fem::quadrature
{
    class Tetrahedron : public IntegrationRule
    {
    public:
        Tetrahedron(int order)
            : IntegrationRule(order)
        {
        }
        int n_points() const override;
        IntegrationPoint point(int i) const override;
    };
}