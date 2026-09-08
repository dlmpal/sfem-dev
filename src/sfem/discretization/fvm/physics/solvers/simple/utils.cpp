#include "utils.hpp"
#include <sfem/discretization/fvm/physics/kernels/transient.hpp>
#include <sfem/discretization/fvm/physics/kernels/convection.hpp>
#include <sfem/discretization/fvm/physics/kernels/laplacian.hpp>
#include <sfem/discretization/fvm/physics/kernels/source.hpp>
#include <sfem/mesh/utils/loop_utils.hpp>

namespace sfem::fvm::simple
{
    //=============================================================================
    void setup_momentum_eqn(FVField u, FVField P_, Field &mu,
                            std::vector<real_t> &mdot,
                            Equation &eqn, int dir)
    {
        eqn.add_kernel(Convection(u, mdot));
        eqn.add_kernel(Laplacian(u, mu));
        auto rhs = [P_, dir](const FVField &, int cell_idx, std::span<real_t> dpdxi)
        {
            dpdxi[0] = -P_.cell_grad(cell_idx)(dir);
        };
        eqn.add_kernel(Source(u, rhs));
    }
    //=============================================================================
    void setup_pressure_eqn(FVField &P, Field &D, Field &rho,
                            const std::vector<real_t> &mdot,
                            Equation &eqn)
    {
        eqn.add_kernel(Laplacian(P, D));

        auto pressure_rhs = [&](la::MatSet, la::VecSet b)
        {
            const auto V = P.space();
            auto work = [&](const mesh::Mesh &,
                            const mesh::Region &,
                            const mesh::Cell &,
                            int facet_idx)
            {
                const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
                const auto [owner, neighbour] = adjacent_cells;
                std::array<real_t, 2> values = {};
                if (owner == neighbour)
                {
                    values[0] = -mdot[facet_idx] / rho.cell_value(facet_idx);
                }
                else
                {
                    values[0] = -mdot[facet_idx] / rho.facet_value(facet_idx);
                    values[1] = mdot[facet_idx] / rho.facet_value(facet_idx);
                }
                b(adjacent_cells, values);
            };
            mesh::utils::for_all_facets(*V->mesh(), work);
        };
        eqn.add_kernel(pressure_rhs);
    }
    //=============================================================================
    void compute_mass_flux(const std::vector<FVField> &U,
                           const FVField &P,
                           const FVField &D,
                           const Field &rho,
                           std::vector<real_t> &mdot)
    {
        const auto V = P.space();
        const auto mesh = V->mesh();
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &region,
                        const mesh::Cell &,
                        int facet_idx)
        {
            const auto [owner, neighbour] = V->facet_adjacent_cells(facet_idx);
            const real_t g = V->facet_interp_factor(facet_idx);
            const geo::Vec3 Sf = V->facet_area_vec(facet_idx);

            mdot[facet_idx] = 0.0;
            if (owner != neighbour)
            {
                const real_t rhof = rho.facet_value(facet_idx);
                for (int dir = 0; dir < mesh->pdim(); dir++)
                {
                    mdot[facet_idx] += rhof * U[dir].facet_value(facet_idx) * Sf(dir);
                }

                // Rhie-Chow correction
                const real_t Df = D.facet_value(facet_idx);
                const geo::Vec3 gradPf = P.facet_grad(facet_idx);
                const geo::Vec3 gradPf_interp = g * P.cell_grad(owner) + (1 - g) * P.cell_grad(neighbour);
                mdot[facet_idx] += -rhof * Df * geo::inner(gradPf - gradPf_interp, Sf);
            }
            else
            {
                /// @todo Handle other BC types?
                /// Can we just use U_[dir].facet_value(facet_idx)?
                for (int dir = 0; dir < mesh->pdim(); dir++)
                {
                    const FVBC &ubc = U[dir].boundary_condition();
                    real_t uf = U[dir].cell_value(owner);
                    const real_t rhof = rho.cell_value(owner);
                    if (ubc.region_type(region.name()) == fvm::BCType::dirichlet)
                    {
                        uf = ubc.value(facet_idx);
                    }
                    mdot[facet_idx] += rhof * uf * Sf(dir);
                }
            }
        };
        mesh::utils::for_all_facets(*mesh, work);
    }
    //=============================================================================
    real_t compute_mass_residual(std::shared_ptr<const FVSpace> V,
                                 const std::vector<real_t> &flux_)
    {
        // Quick access
        const auto mesh = V->mesh();

        // Net mass flux out of every cell
        la::Vector imbalance(V->index_map(), 1);
        auto add_values = la::create_vecset(imbalance);

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int facet_idx)
        {
            const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
            const auto [owner, neighbour] = adjacent_cells;
            if (owner == neighbour)
            {
                const std::array<int, 1> idx = {owner};
                const std::array<real_t, 1> values = {flux_[facet_idx]};
                add_values(idx, values);
            }
            else
            {
                const std::array<real_t, 2> values = {flux_[facet_idx],
                                                      -flux_[facet_idx]};
                add_values(adjacent_cells, values);
            }
        };
        mesh::utils::for_all_facets(*mesh, work);
        imbalance.assemble();

        return la::norm(imbalance, la::NormType::l1);;
    }
}