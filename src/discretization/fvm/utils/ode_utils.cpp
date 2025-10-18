#include "ode_utils.hpp"

namespace sfem::fvm::ode
{
    //=============================================================================
    RHSFunction create_rhs(std::shared_ptr<const fvm::FVField> phi,
                           std::shared_ptr<const fvm::NumericalFlux> nflux,
                           SourceFunction src)
    {
        return [=](const la::Vector &S,
                   la::Vector &rhs,
                   real_t time)
        {
            // Quick access
            const auto V = phi->space();
            const auto mesh = V->mesh();
            const auto topology = mesh->topology();
            const auto flux = nflux->flux_function();

            // Reset RHS vector
            rhs.set_all(0.0);

            // Flux contribution
            std::vector<real_t> state_left(flux->n_comp());
            std::vector<real_t> state_right(flux->n_comp());
            std::vector<real_t> normal_flux(flux->n_comp());

            for (const auto &region : mesh->regions())
            {
                // Skip boundary regions
                if (region.dim() < mesh->pdim())
                {
                    continue;
                }

                for (const auto &[facet, facet_idx] : mesh->region_facets(region.name()))
                {
                    // Only integrate locally owned and non-boundary facets
                    if (topology->entity_index_map(mesh->pdim() - 1)->is_ghost(facet_idx) ||
                        V->is_boundary(facet_idx))
                    {
                        continue;
                    }

                    // Get the facet's adjacent cells
                    const auto [cell_idx_left, cell_idx_right] = V->facet_adjacent_cells(facet_idx);

                    // Get the adjacent cell states
                    for (int i = 0; i < flux->n_comp(); i++)
                    {
                        state_left[i] = S(cell_idx_left, i);
                        state_right[i] = S(cell_idx_right, i);
                    }

                    // Facet area and normal vector
                    const real_t area = V->facet_area(facet_idx);
                    const auto normal = V->facet_normal(facet_idx);

                    // Compute the (numerical) normal flux at the face
                    nflux->compute_normal_flux(state_left, state_right, normal, normal_flux);

                    // Inverse left and right cell volumes
                    const real_t vol_inv_left = 1 / V->cell_volume(cell_idx_left);
                    const real_t vol_inv_right = 1 / V->cell_volume(cell_idx_right);

                    // Add the flux contributions to the RHS vector
                    for (int i = 0; i < flux->n_comp(); i++)
                    {
                        rhs(cell_idx_left, i) -= normal_flux[i] * area * vol_inv_left;
                        rhs(cell_idx_right, i) -= -normal_flux[i] * area * vol_inv_right;
                    }
                }

                // Source term contribution
                if (src)
                {
                    std::vector<real_t> source(flux->n_comp());
                    for (const auto &[cell, cell_idx] : mesh->region_cells(region.name()))
                    {
                        // Only integrate locally owned cells
                        if (topology->entity_index_map(mesh->pdim())->is_ghost(cell_idx))
                        {
                            continue;
                        }

                        src(V->cell_midpoint(cell_idx), source, time);
                        for (int i = 0; i < flux->n_comp(); i++)
                        {
                            rhs(cell_idx, i) += source[i];
                        }
                    }
                }
            }

            rhs.assemble();
        };
    }
}