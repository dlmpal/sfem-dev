#include "fv_function.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FVFunction::FVFunction(std::shared_ptr<const FVSpace> fv_space,
                           const std::vector<std::string> &components)
        : Function(fv_space->index_map(), components),
          fv_space_(fv_space)
    {
    }
    //=============================================================================
    std::shared_ptr<const FVSpace> FVFunction::space() const
    {
        return fv_space_;
    }
}