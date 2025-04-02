#pragma once

#include "quadrature.hpp"

namespace sfem::fem::quadrature
{
    class Gauss1D : public IntegrationRule
    {
    public:
        Gauss1D(int order)
            : IntegrationRule(order)
        {
        }
        int n_points() const override;
        IntegrationPoint point(int i) const override;
    };

    class Gauss2D : public IntegrationRule
    {
    public:
        Gauss2D(int order)
            : IntegrationRule(order)
        {
        }
        int n_points() const override;
        IntegrationPoint point(int i) const override;
    };

    class Gauss3D : public IntegrationRule
    {
    public:
        Gauss3D(int order)
            : IntegrationRule(order)
        {
        }
        int n_points() const override;
        IntegrationPoint point(int i) const override;
    };
}