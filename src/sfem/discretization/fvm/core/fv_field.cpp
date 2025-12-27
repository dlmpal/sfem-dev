#include "fv_field.hpp"
#include <sfem/discretization/fvm/core/fv_bc.hpp>
#include <sfem/discretization/fvm/core/fv_gradient.hpp>
#include <sfem/la/native/vector.hpp>

namespace sfem::fvm
{
    //=============================================================================
    IField::IField(const std::vector<std::string> &components)
        : components_(components)
    {
    }
    //=============================================================================
    std::vector<std::string> IField::components() const
    {
        return components_;
    }
    //=============================================================================
    int IField::n_comp() const
    {
        return static_cast<int>(components_.size());
    }
    //=============================================================================
    ConstantField::ConstantField(const std::vector<std::string> &components,
                                 const std::vector<real_t> &value)
        : IField(components),
          value_(value)
    {
    }
    //=============================================================================
    ConstantField::ConstantField(const std::string &name, real_t value)
        : ConstantField(std::vector<std::string>{name},
                        std::vector<real_t>{value})
    {
    }
    //=============================================================================
    real_t &ConstantField::cell_value(int, int comp_idx)
    {
        return value_[comp_idx];
    }
    //=============================================================================
    real_t ConstantField::cell_value(int, int comp_idx) const
    {
        return value_[comp_idx];
    }
    //=============================================================================
    real_t ConstantField::facet_value(int, int comp_idx) const
    {
        return value_[comp_idx];
    }
    //=============================================================================
    geo::Vec3 ConstantField::cell_grad(int, int) const
    {
        return {};
    }
    //=============================================================================
    geo::Vec3 ConstantField::facet_grad(int, int) const
    {
        return {};
    }
    //=============================================================================
    FVField::FVField(std::shared_ptr<const FVSpace> V,
                     const std::vector<std::string> &components,
                     GradientMethod gradient_method)
        : IField(components),
          V_(V),
          topo_(V_->mesh()->topology()),
          bc_(std::make_shared<FVBC>(*V, n_comp())),
          values_(std::make_shared<la::Vector>(V->index_map(), n_comp())),
          gradient_method_(gradient_method)
    {
        if (gradient_method_ != GradientMethod::none)
        {
            const int dim = V->mesh()->pdim();
            grad_ = std::make_shared<la::Vector>(V->index_map(), n_comp() * dim);
        }
    }
    //=============================================================================
    std::shared_ptr<const FVSpace> FVField::space() const
    {
        return V_;
    }
    //=============================================================================
    FVBC &FVField::boundary_condition()
    {
        return *bc_;
    }
    //=============================================================================
    const FVBC &FVField::boundary_condition() const
    {
        return *bc_;
    }
    //=============================================================================
    la::Vector &FVField::values()
    {
        return *values_;
    }
    //=============================================================================
    const la::Vector &FVField::values() const
    {
        return *values_;
    }
    //=============================================================================
    GradientMethod FVField::grad_method() const
    {
        return gradient_method_;
    }
    //=============================================================================
    la::Vector &FVField::grad()
    {
        return *grad_;
    }
    //=============================================================================
    const la::Vector &FVField::grad() const
    {
        return *grad_;
    }
    //=============================================================================
    real_t &FVField::cell_value(int cell_idx, int comp_idx)
    {
        return (*values_)(cell_idx, comp_idx);
    }
    //=============================================================================
    real_t FVField::cell_value(int cell_idx, int comp_idx) const
    {
        return (*values_)(cell_idx, comp_idx);
    }
    //=============================================================================
    real_t FVField::facet_value(int facet_idx, int comp_idx) const
    {
        const auto [owner, neighbour] = V_->facet_adjacent_cells(facet_idx);
        if (owner == neighbour)
        {
            const int tag = topo_->facets()[facet_idx].tag;
            const auto region = V_->mesh()->get_region_by_tag(tag);
            const auto bc_type = bc_->region_type(region.name());

            const real_t phiP = cell_value(owner);
            const real_t dPf = V_->facet_cell_distances(facet_idx)[0];

            if (bc_type == BCType::dirichlet)
            {
                return bc_->value(facet_idx, comp_idx);
            }
            else if (bc_type == BCType::neumann)
            {
                return phiP - dPf * bc_->value(facet_idx, comp_idx);
            }
            else if (bc_type == BCType::robin)
            {
            }
            else // Zero Neumann
            {
                return phiP;
            }
        }
        else
        {
            const real_t g = V_->facet_interp_factor(facet_idx);
            const real_t phiP = cell_value(owner, comp_idx);
            const real_t phiN = cell_value(neighbour, comp_idx);
            return g * phiP + (1 - g) * phiN;
        }
    }
    //=============================================================================
    geo::Vec3 FVField::cell_grad(int cell_idx, int comp_idx) const
    {
        if (gradient_method_ == GradientMethod::none)
        {
            return geo::Vec3();
        }
        else
        {
            const int dim = V_->mesh()->pdim();
            geo::Vec3 grad;
            for (int dir = 0; dir < dim; dir++)
            {
                grad(dir) = (*grad_)(cell_idx, comp_idx * dim + dir);
            }
            return grad;
        }
    }
    //=============================================================================
    geo::Vec3 FVField::facet_grad(int facet_idx, int comp_idx) const
    {
        if (gradient_method_ == GradientMethod::none)
        {
            return geo::Vec3();
        }

        const auto [owner, neighbour] = V_->facet_adjacent_cells(facet_idx);
        if (owner == neighbour)
        {
            /// @todo Use BC info
            return cell_grad(owner, comp_idx);
        }
        else
        {
            const geo::Vec3 dPN = V_->facet_intercell_distance(facet_idx);
            const geo::Vec3 ePN = dPN.normalize();
            const real_t g = V_->facet_interp_factor(facet_idx);
            const real_t phiP = cell_value(owner, comp_idx);
            const real_t phiN = cell_value(neighbour, comp_idx);
            const geo::Vec3 gradP = cell_grad(owner, comp_idx);
            const geo::Vec3 gradN = cell_grad(neighbour, comp_idx);
            const geo::Vec3 grad_avg = g * gradP + (1 - g) * gradN;
            return grad_avg + ePN * ((phiN - phiP) / dPN.mag() - geo::inner(grad_avg, ePN));
        }
    }
    //=============================================================================
    void FVField::update_gradient()
    {
        if (gradient_method_ == GradientMethod::none)
        {
            return;
        }
        else if (gradient_method_ == GradientMethod::green_gauss)
        {
            green_gauss_gradient(*this);
        }
        else if (gradient_method_ == GradientMethod::least_squares)
        {
            least_squares_gradient(*this);
        }
    }
}