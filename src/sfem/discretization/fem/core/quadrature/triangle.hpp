#pragma once

#include <sfem/discretization/fem/core/quadrature/quadrature.hpp>

namespace sfem::fem::quadrature
{
    class Triangle : public IntegrationRule
    {
    public:
        Triangle(int n_points);
        real_t weight(int qpt_idx) const override;
        std::array<real_t, 3> point(int qpt_idx) const override;
    };
}