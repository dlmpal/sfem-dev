#pragma once

#include "quadrature.hpp"

namespace sfem::fem::quadrature
{
    class Triangle : public IntegrationRule
    {
    public:
        Triangle(int n_points);
        IntegrationPoint point(int i) const override;
    };
}