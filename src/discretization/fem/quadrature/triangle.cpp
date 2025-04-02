#include "triangle.hpp"

namespace sfem::fem::quadrature
{
    //=============================================================================
    int Triangle::n_points() const
    {
        return 0;
    }
    //=============================================================================
    IntegrationPoint Triangle::point(int i) const
    {
        return IntegrationPoint{};
    }
}