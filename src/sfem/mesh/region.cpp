#include "region.hpp"
#include <sfem/base/error.hpp>

namespace sfem::mesh
{
    //=============================================================================
    Region::Region(const std::string &name, int tag, int dim)
        : name_(name),
          tag_(tag),
          dim_(dim)
    {
        if (tag_ < 0)
        {
            SFEM_ERROR(std::format("Region {} has invalid tag {} (<0)\n", name_, dim_));
        }
        if (dim_ < 0)
        {
            SFEM_ERROR(std::format("Region {} has invalid dimension of {} (<0)\n", name_, dim_));
        }
    }
    //=============================================================================
    std::string Region::name() const
    {
        return name_;
    }
    //=============================================================================
    int Region::tag() const
    {
        return tag_;
    }
    //=============================================================================
    int Region::dim() const
    {
        return dim_;
    }
}
