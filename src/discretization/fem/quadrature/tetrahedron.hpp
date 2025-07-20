#pragma once

#include "quadrature.hpp"

namespace sfem::fem::quadrature
{
    class Tetrahedron : public IntegrationRule
    {
    public:
        Tetrahedron(int n_points);
        IntegrationPoint point(int i) const override;
    };
}