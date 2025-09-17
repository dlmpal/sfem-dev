#include "function.hpp"

namespace sfem
{
  //=============================================================================
  Function::Function(std::shared_ptr<const IndexMap> index_map,
                     const std::vector<std::string> &components)
      : la::Vector(index_map, static_cast<int>(components.size())),
        components_(components)
  {
  }
  //=============================================================================
  std::vector<std::string> Function::components() const
  {
    return components_;
  }
  //=============================================================================
  int Function::n_comp() const
  {
    return bs_;
  }
  //=============================================================================
  int Function::comp_idx(const std::string &component) const
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
      log_msg(msg, LogLevel::warning);
      return -1;
    }
    else
    {
      return static_cast<int>(std::distance(components_.cbegin(), it));
    }
  }
  //=============================================================================
  MeshFunction::MeshFunction(std::shared_ptr<const mesh::Mesh> mesh, int dim,
                             const std::vector<std::string> &components)
      : Function(mesh->topology()->entity_index_map(dim), components),
        mesh_(mesh)
  {
  }
  //=============================================================================
  std::shared_ptr<const mesh::Mesh> MeshFunction::mesh() const
  {
    return mesh_;
  }
}