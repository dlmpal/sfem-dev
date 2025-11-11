#include "convection.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void _convection_impl(const FVField &phi,
                          const FVBC &bc,
                          const Coefficient &vel,
                          la::MatSet lhs, la::VecSet rhs)
    {
        // Finite volume space
        const auto V = phi.space();
        const int dim = V->mesh()->pdim();

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &region,
                        const mesh::Cell &,
                        int facet_idx)
        {
            // Cells adjacent to the facet
            const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
            const auto &[cell_idx1, cell_idx2] = adjacent_cells;

            // Facet area, interpolation factor and normal vector
            const real_t area = V->facet_area(facet_idx);
            const real_t g = V->facet_interp_factor(facet_idx);
            const geo::Vec3 normal = V->facet_normal(facet_idx);

            // Compute the normal velocity
            real_t normal_vel = 0.0;
            for (int i = 0; i < dim; i++)
            {
                normal_vel += (g * vel(cell_idx1, i) + (1 - g) * vel(cell_idx2, i)) * normal(i);
            }

            // Compute the flux
            const real_t flux = normal_vel * area;

            // Boundary facets
            if (cell_idx1 == cell_idx2)
            {
                const std::array<int, 1> idx = {cell_idx1};
                std::array<real_t, 1> lhs_value{};
                std::array<real_t, 1> rhs_value{};

                if (bc.types_.at(region.name()) == BCType::dirichlet)
                {
                    if (normal_vel > 0)
                    {
                        lhs_value[0] = flux;
                    }
                    else
                    {
                        const real_t phi_facet = bc.values_.at(facet_idx)[0];
                        rhs_value[0] = -flux * phi_facet;
                    }
                }
                else if (bc.types_.at(region.name()) == BCType::neumann)
                {
                    if (normal_vel > 0)
                    {
                        lhs_value[0] = flux;
                    }
                    else
                    {
                        const real_t grad_facet = bc.values_.at(facet_idx)[0];
                        const real_t d = V->facet_cell_distances(facet_idx)[0];
                        lhs_value[0] = flux;
                        rhs_value[0] = -flux * d * grad_facet;
                    }
                }
                else if (bc.types_.at(region.name()) == BCType::zero_neumann)
                {
                    lhs_value[0] = flux;
                }

                lhs(idx, idx, lhs_value);
                rhs(idx, rhs_value);
            }
            // Internal facets
            else
            {
                // Upwind differencing
                const real_t w = normal_vel > 0 ? 1.0 : 0.0;

                const std::array<int, 2> idxs = {cell_idx1, cell_idx2};
                const std::array<real_t, 4> lhs_values = {w * flux, (1 - w) * flux,
                                                          -w * flux, -(1 - w) * flux};
                lhs(idxs, idxs, lhs_values);
            }
        };
        mesh::utils::for_all_facets(*V->mesh(), work);
    }
    //=============================================================================
    void convection(const FVField &phi,
                    const FVBC &bc,
                    const Coefficient &vel,
                    la::MatSet lhs, la::VecSet rhs)
    {
        SFEM_CHECK_SIZES(phi.space()->mesh()->pdim(), vel.n_comp());

        _convection_impl(phi, bc, vel, lhs, rhs);

        // switch (scheme)
        // {
        // case DifferencingScheme::upwind:
        //     _convection_impl(phi, grad, bc, vel, lhs, rhs, upwind, implicit);
        //     break;
        // case DifferencingScheme::linear:
        //     _convection_impl(phi, grad, bc, vel, lhs, rhs, linear, implicit);
        //     break;
        // default:
        //     SFEM_ERROR(std::format("Invalid differencing scheme: {}\n",
        //                            static_cast<int>(scheme)));
        //     break;
        // }
    }
}