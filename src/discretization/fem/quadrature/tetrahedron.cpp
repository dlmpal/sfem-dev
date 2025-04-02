#include "tetrahedron.hpp"

namespace sfem::fem::quadrature
{
    //=============================================================================
    int Tetrahedron::n_points() const
    {
        return 0;
    }
    //=============================================================================
    IntegrationPoint Tetrahedron::point(int i) const
    {
        return IntegrationPoint{};
    }
}