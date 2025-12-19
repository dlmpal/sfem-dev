#include "quadrature.hpp"

namespace sfem::fem
{
    //=============================================================================
    IntegrationRule::IntegrationRule(int n_points)
        : n_points_(n_points)
    {
    }
    //=============================================================================
    int IntegrationRule::n_points() const
    {
        return n_points_;
    }
    //=============================================================================
    void IntegrationRule::set_n_points(int n_points)
    {
        n_points_ = n_points;
    }
}