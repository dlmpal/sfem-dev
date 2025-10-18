#include "field.hpp"
#include "../base/logging.hpp"
#include <format>

namespace sfem
{
  //=============================================================================
  Field::Field(std::shared_ptr<const IndexMap> index_map,
               const std::vector<std::string> &components)
      : la::Vector(index_map, static_cast<int>(components.size())),
        components_(components)
  {
  }
  //=============================================================================
  std::vector<std::string> Field::components() const
  {
    return components_;
  }
  //=============================================================================
  int Field::n_comp() const
  {
    return bs_;
  }
  //=============================================================================
  int Field::comp_idx(const std::string &component) const
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
}