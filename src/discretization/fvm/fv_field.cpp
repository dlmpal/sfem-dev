#include "fv_field.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FVField::FVField(std::shared_ptr<const FVSpace> fv_space,
                     const std::vector<std::string> &components)
        : Field(fv_space->index_map(), components),
          fv_space_(fv_space)
    {
    }
    //=============================================================================
    std::shared_ptr<const FVSpace> FVField::space() const
    {
        return fv_space_;
    }
}