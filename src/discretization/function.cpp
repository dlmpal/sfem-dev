#include "function.hpp"

namespace sfem
{
    //=============================================================================
    Function::Function(int n_comp)
        : n_comp_(n_comp)
    {
    }
    //=============================================================================
    int Function::n_comp() const
    {
        return n_comp_;
    }
    //=============================================================================
    ConstantFunction::ConstantFunction(real_t value)
        : Function(1), value_(value)
    {
    }
    //=============================================================================
    real_t &ConstantFunction::value()
    {
        return value_;
    }
    //=============================================================================
    real_t ConstantFunction::value() const
    {
        return value_;
    }
    //=============================================================================
    real_t &ConstantFunction::operator()([[maybe_unused]] int idx, [[maybe_unused]] int comp)
    {
        return value_;
    }
    //=============================================================================

    real_t ConstantFunction::operator()([[maybe_unused]] int idx, [[maybe_unused]] int comp) const
    {
        return value_;
    }
    //=============================================================================
    MeshFunction::MeshFunction(std::shared_ptr<mesh::Mesh> mesh, int dim, int n_comp)
        : Function(n_comp),
          mesh_(mesh),
          values_(mesh_->topology()->entity_index_map(dim), n_comp)
    {
    }
    //=============================================================================
    std::shared_ptr<const mesh::Mesh> MeshFunction::mesh() const
    {
        return mesh_;
    }
    //=============================================================================
    la::Vector &MeshFunction::values()
    {
        return values_;
    }
    //=============================================================================
    const la::Vector &MeshFunction::values() const
    {
        return values_;
    }
    //=============================================================================
    real_t &MeshFunction::operator()(int idx, int comp)
    {
        return values_(idx, comp);
    }
    //=============================================================================
    real_t MeshFunction::operator()(int idx, int comp) const
    {
        return values_(idx, comp);
    }
}