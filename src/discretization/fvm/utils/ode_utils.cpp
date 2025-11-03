#include "ode_utils.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm::ode
{
    //=============================================================================
    RHSFunction create_rhs(std::shared_ptr<const fvm::FVField> phi,
                           std::shared_ptr<const fvm::NumericalFlux> nflux,
                           FieldFunction src)
    {
        return [=](const la::Vector &S,
                   la::Vector &rhs,
                   real_t time)
        {
            // Reset RHS vector
            rhs.set_all(0.0);

            // Finite volume space and flux function
            const auto V = phi->space();
            const auto flux = nflux->flux_function();

            // Left and right state fluxes and normal flux
            std::vector<real_t> state1(flux->n_comp());
            std::vector<real_t> state2(flux->n_comp());
            std::vector<real_t> normal_flux(flux->n_comp());

            // Add flux contribution
            auto facet_work = [&](const mesh::Mesh &,
                                  const mesh::Region &,
                                  const mesh::Cell &,
                                  int facet_idx)
            {
                // Get the facet's adjacent cells
                const auto [cell_idx1, cell_idx2] = V->facet_adjacent_cells(facet_idx);

                // Get the adjacent cell states
                for (int i = 0; i < flux->n_comp(); i++)
                {
                    state1[i] = S(cell_idx1, i);
                    state2[i] = S(cell_idx2, i);
                }

                // Facet area and normal vector
                const real_t area = V->facet_area(facet_idx);
                const auto normal = V->facet_normal(facet_idx);

                // Compute the (numerical) normal flux at the face
                nflux->compute_normal_flux(state1, state2, normal, normal_flux);

                // Inverse left and right cell volumes
                const real_t vol_inv1 = 1 / V->cell_volume(cell_idx1);
                const real_t vol_inv2 = 1 / V->cell_volume(cell_idx2);

                // Add the flux contributions to the RHS vector
                for (int i = 0; i < flux->n_comp(); i++)
                {
                    rhs(cell_idx1, i) -= normal_flux[i] * area * vol_inv1;

                    /// @todo cleanup
                    if (cell_idx1 != cell_idx2)
                    {
                        rhs(cell_idx2, i) -= -normal_flux[i] * area * vol_inv2;
                    }
                }
            };
            mesh::utils::for_all_facets(*V->mesh(), facet_work);

            // Add source term contribution
            if (src)
            {
                std::vector<real_t> source(flux->n_comp());
                auto cell_work = [&](const mesh::Mesh &,
                                     const mesh::Region &,
                                     const mesh::Cell &,
                                     int cell_idx)
                {
                    src(V->cell_midpoint(cell_idx), source, time);
                    for (int i = 0; i < flux->n_comp(); i++)
                    {
                        rhs(cell_idx, i) += source[i];
                    }
                };
                mesh::utils::for_all_cells(*V->mesh(), cell_work);
            };

            // RHS was incrementally constructed - needs assembly
            rhs.assemble();
        };
    }
}