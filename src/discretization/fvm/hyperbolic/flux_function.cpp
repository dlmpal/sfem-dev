#include "flux_function.hpp"
#include <cmath>

namespace sfem::fvm
{
    //=============================================================================
    FluxFunction::FluxFunction(int n_comp, int dim)
        : n_comp_(n_comp),
          dim_(dim)
    {
    }
    //=============================================================================
    int FluxFunction::n_comp() const
    {
        return n_comp_;
    }
    //=============================================================================
    int FluxFunction::dim() const
    {
        return dim_;
    }
}