#include "fe_field.hpp"

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
    real_t ConstantField::value(int, mesh::CellType,
                                const std::array<real_t, 3> &, int comp_idx) const
    {
        return value_[comp_idx];
    }
    //=============================================================================
    geo::Vec3 ConstantField::grad(int, mesh::CellType,
                                  const std::array<real_t, 3> &, int) const
    {
        return geo::Vec3{0, 0, 0};
    }
    //=============================================================================
    FEField::FEField(std::shared_ptr<const FESpace> V,
                     const std::vector<std::string> &components)
        : Field(components),
          V_(V),
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
    real_t FEField::value(int cell_idx, mesh::CellType cell_type,
                          const std::array<real_t, 3> &pt, int comp_idx) const
    {
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
    geo::Vec3 FEField::grad(int cell_idx, mesh::CellType cell_type,
                            const std::array<real_t, 3> &pt, int comp_idx) const
    {
        const auto element = V_->element(cell_type);
        const auto elem_dof = V_->cell_dof(cell_idx);
        const auto elem_pts = V_->cell_dof_points(cell_idx);
        const auto data = element->transform(element->dim(), pt, elem_pts);

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
}