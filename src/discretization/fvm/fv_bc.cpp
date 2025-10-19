#include "fv_bc.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FVBC::FVBC(std::shared_ptr<const FVField> phi)
        : phi_(phi)
    {
        // Set all boundary region BC types to
        // zero Neumann by default
        const auto mesh = phi_->space()->mesh();
        for (const auto &region : mesh->regions())
        {
            if (region.dim() < mesh->pdim())
            {
                types_[region.name()] = BCType::zero_neumann;
            }
        }
    }
    //=============================================================================
    void FVBC::set_value(const std::string &region,
                         const std::string &comp_name,
                         BCType type, std::array<real_t, 2> values)
    {
        // Quick access
        const auto V = phi_->space();
        const auto mesh = V->mesh();
        const int n_comp = phi_->n_comp();
        const int comp_idx = phi_->comp_idx(comp_name);

        // Set region BC type
        types_[region] = type;

        // Loop over region facets and store BC value for boundary cells
        for (const auto &[facet, facet_idx] : mesh->region_facets(region))
        {
            values_[facet_idx * n_comp + comp_idx][0] = values[0];
            values_[facet_idx * n_comp + comp_idx][1] = values[1];
        }
    }
}