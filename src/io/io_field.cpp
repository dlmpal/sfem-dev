#include "io_field.hpp"
#include "../base/logging.hpp"
#include <format>

namespace sfem::io
{
    //=============================================================================
    IOField::IOField(const std::vector<std::string> &components, ValueAccessor func)
        : components_(components),
          values_(func)
    {
    }
    //=============================================================================
    std::vector<std::string> IOField::components() const
    {
        return components_;
    }
    //=============================================================================
    int IOField::comp_idx(const std::string &component) const
    {
        auto it = std::find(components_.cbegin(),
                            components_.end(),
                            component);
        if (it == components_.end())
        {
            auto msg = std::format("Component {} not found in: [ ", component);
            for (const auto &comp : components_)
            {
                msg += std::format("{} ", comp);
            }
            msg += "]\n";
            log_msg(msg, false, LogLevel::warning);
            return -1;
        }
        else
        {
            return static_cast<int>(std::distance(components_.cbegin(), it));
        }
    }
    //=============================================================================
    real_t IOField::operator()(int idx, int comp_idx) const
    {
        return values_(idx, comp_idx);
    }
}