#include "ode_utils.hpp"
#include <sfem/la/native/vector.hpp>
#include <sfem/mesh/utils/loop_utils.hpp>

namespace sfem::fvm::ode
{
    //=============================================================================
    RHSFunction create_rhs(std::shared_ptr<const fvm::FVField> phi,
                           std::shared_ptr<const fvm::NumericalFlux> nflux)
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
            std::vector<real_t> uP(flux->n_comp());
            std::vector<real_t> uN(flux->n_comp());
            std::vector<real_t> normal_flux(flux->n_comp());

            // Add flux contribution
            auto facet_work = [&](const mesh::Mesh &,
                                  const mesh::Region &,
                                  const mesh::Cell &,
                                  int facet_idx)
            {
                // Get the facet's adjacent cells
                const auto [owner, neighbour] = V->facet_adjacent_cells(facet_idx);

                // Get the adjacent cell states
                for (int i = 0; i < flux->n_comp(); i++)
                {
                    uP[i] = S(owner, i);
                    uN[i] = S(neighbour, i);
                }

                // Facet area vector and area
                const geo::Vec3 Sf = V->facet_area_vec(facet_idx);
                const real_t Af = Sf.mag();

                // Compute the (numerical) normal flux at the face
                nflux->compute_normal_flux(uP, uN, Sf.normalize(), normal_flux);

                // Inverse left and right cell volumes
                const real_t vol_inv1 = 1 / V->cell_volume(owner);
                const real_t vol_inv2 = 1 / V->cell_volume(neighbour);

                // Add the flux contributions to the RHS vector
                for (int i = 0; i < flux->n_comp(); i++)
                {
                    rhs(owner, i) -= normal_flux[i] * Af * vol_inv1;

                    /// @todo cleanup
                    if (owner != neighbour)
                    {
                        rhs(neighbour, i) -= -normal_flux[i] * Af * vol_inv2;
                    }
                }
            };
            mesh::utils::for_all_facets(*V->mesh(), facet_work);

            // RHS was incrementally constructed - needs assembly
            rhs.assemble();
        };
    }
}