#include "fe_field.hpp"
#include "../../mesh/utils/geo_utils.hpp"

namespace sfem::fem
{
    //=============================================================================
    Field::Field(const std::vector<std::string> &components)
        : components_(components)
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
        return static_cast<int>(components_.size());
    }
    //=============================================================================
    ConstantField::ConstantField(const std::vector<std::string> &components,
                                 const std::vector<real_t> &value)
        : Field(components),
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
    real_t ConstantField::cell_value(int, const std::array<real_t, 3> &, int comp_idx) const
    {
        return value_[comp_idx];
    }
    //=============================================================================
    real_t ConstantField::facet_value(int, const std::array<real_t, 3> &, int comp_idx) const
    {
        return value_[comp_idx];
    }
    //=============================================================================
    geo::Vec3 ConstantField::cell_grad(int, const std::array<real_t, 3> &, int) const
    {
        return geo::Vec3{0, 0, 0};
    }
    //=============================================================================
    geo::Vec3 ConstantField::facet_grad(int, const std::array<real_t, 3> &, int) const
    {
        return geo::Vec3{0, 0, 0};
    }
    //=============================================================================
    FEField::FEField(std::shared_ptr<const FESpace> V,
                     const std::vector<std::string> &components)
        : Field(components),
          V_(V),
          topo_(V_->mesh()->topology()),
          dof_values_(std::make_shared<la::Vector>(V->index_map(), n_comp()))
    {
    }
    //=============================================================================
    std::shared_ptr<const FESpace> FEField::space() const
    {
        return V_;
    }
    //=============================================================================
    la::Vector &FEField::dof_values()
    {
        return *dof_values_;
    }
    //=============================================================================
    const la::Vector &FEField::dof_values() const
    {
        return *dof_values_;
    }
    //=============================================================================
    void FEField::cell_values(int cell_idx, std::span<real_t> values) const
    {
        const auto elem_dof = V_->cell_dof(cell_idx);
        const int n_nodes = static_cast<int>(elem_dof.size());
        const int n_comp_ = n_comp();
        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < n_comp_; j++)
            {
                values[i * n_comp_ + j] = (*dof_values_)(elem_dof[i], j);
            }
        }
    }
    //=============================================================================
    void FEField::cell_values(int cell_idx, int comp_idx,
                              std::span<real_t> values) const
    {
        const auto elem_dof = V_->cell_dof(cell_idx);
        const int n_nodes = static_cast<int>(elem_dof.size());
        for (int i = 0; i < n_nodes; i++)
        {
            values[i] = (*dof_values_)(elem_dof[i], comp_idx);
        }
    }
    //=============================================================================
    real_t FEField::cell_value(int cell_idx, const std::array<real_t, 3> &pt, int comp_idx) const
    {
        const auto cell_type = topo_->cells().at(cell_idx).type;
        const auto element = V_->element(cell_type);
        const auto elem_dof = V_->cell_dof(cell_idx);

        la::DenseMatrix N(element->n_nodes(), 1);
        element->eval_shape(pt, N);

        real_t value = 0.0;
        for (int i = 0; i < element->n_nodes(); i++)
        {
            value += (*dof_values_)(elem_dof[i], comp_idx) * N(i, 0);
        }
        return value;
    }
    //=============================================================================
    real_t FEField::facet_value(int facet_idx, const std::array<real_t, 3> &pt, int comp_idx) const
    {
        const int dim = topo_->dim();
        const int cell_idx = topo_->entity_owner(facet_idx, dim - 1);
        const int rel_idx = topo_->entity_rel_idx(cell_idx, dim, facet_idx, dim - 1);
        const auto cell_type = topo_->cells().at(cell_idx).type;

        const auto pt_cell = mesh::map_facet_to_cell_ref(cell_type, rel_idx, pt);

        return cell_value(cell_idx, pt_cell, comp_idx);
    }
    //=============================================================================
    geo::Vec3 FEField::cell_grad(int cell_idx, const std::array<real_t, 3> &pt, int comp_idx) const
    {
        const auto cell_type = topo_->cells()[cell_idx].type;
        const auto element = V_->element(cell_type);
        const auto elem_dof = V_->cell_dof(cell_idx);
        const auto elem_pts = V_->cell_dof_points(cell_idx);

        const auto data = element->transform(cell_idx, element->dim(), pt, elem_pts);

        geo::Vec3 grad(0.0, 0.0, 0.0);
        for (int i = 0; i < element->n_nodes(); ++i)
        {
            const real_t ui = (*dof_values_)(elem_dof[i], comp_idx);
            for (int dir = 0; dir < element->dim(); dir++)
            {
                grad(dir) += ui * data.dNdX(i, dir);
            }
        }
        return grad;
    }
    //=============================================================================
    geo::Vec3 FEField::facet_grad(int facet_idx, const std::array<real_t, 3> &pt, int comp_idx) const
    {
        const int dim = topo_->dim();
        const int cell_idx = topo_->entity_owner(facet_idx, dim - 1);
        const mesh::CellType cell_type = topo_->cells()[cell_idx].type;
        const int rel_idx = topo_->entity_rel_idx(cell_idx, dim, facet_idx, dim - 1);

        const auto pt_cell = mesh::map_facet_to_cell_ref(cell_type, rel_idx, pt);

        return cell_grad(cell_idx, pt_cell, comp_idx);
    }
}