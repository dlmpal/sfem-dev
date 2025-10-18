#include "dirichlet_bc.hpp"
#include "../../base/error.hpp"
#include <ranges>

namespace sfem::fem
{
    //=============================================================================
    DirichletBC::DirichletBC(std::shared_ptr<const fem::FEField> phi)
        : phi_(phi)
    {
    }
    //=============================================================================
    void DirichletBC::set_value(const std::string &region_name,
                                const std::string &component,
                                real_t value)
    {
        std::vector<real_t> values = {value};
        set_values(region_name, component, values);
    }
    //=============================================================================
    void DirichletBC::set_values(const std::string &region_name,
                                 const std::string &component,
                                 std::span<const real_t> values)
    {
        // Quick access
        const auto fe_space = phi_->space();

        // Get the DoF belonging to the region.
        // If the region is not included, obtain the DoF first
        if (boundary_dof_.contains(region_name) == false)
        {
            auto region_dof = fe_space->boundary_dof(region_name);
            boundary_dof_.insert({region_name, region_dof});
        }
        const auto &dof = boundary_dof_.at(region_name);

        // Store the specified values
        const int n_comp = phi_->n_comp();
        const int comp_idx = phi_->comp_idx(component);
        if (values.size() == 1) ///< Single value provided for all region DoF
        {
            for (std::size_t i = 0; i < dof.size(); i++)
            {
                data_.insert({dof[i] * n_comp + comp_idx, values[0]});
            }
        }
        else ///< Potentially different values for each region DoF
        {
            SFEM_CHECK_SIZES(dof.size(), values.size());
            for (std::size_t i = 0; i < dof.size(); i++)
            {
                data_.insert({dof[i] * n_comp + comp_idx, values[i]});
            }
        }
    }
    //=============================================================================
    std::pair<std::vector<int>, std::vector<real_t>>
    DirichletBC::get_dofs_values() const
    {
        std::vector<int> idxs(data_.size());      ///< DoF indices
        std::vector<real_t> values(data_.size()); ///< DoF values
        for (auto [i, kv] : std::views::enumerate(data_))
        {
            idxs[i] = kv.first;
            values[i] = kv.second;
        }
        return {idxs, values};
    }
    //=============================================================================
    void DirichletBC::reset_values()
    {
        data_.clear();
    }
}