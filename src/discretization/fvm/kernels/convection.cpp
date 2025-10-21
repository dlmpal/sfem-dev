#include "convection.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    template <typename F>
    concept Interp = requires(F f, real_t flux, real_t w,
                              real_t phi1, real_t phi2) {
        { f(flux, w, phi1, phi2) } -> std::same_as<real_t>;
    };
    //=============================================================================
    real_t upwind(real_t flux, real_t w,
                  real_t phi1, real_t phi2)
    {
        if (flux > 0)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    //=============================================================================
    real_t linear(real_t flux, real_t w,
                  real_t phi1, real_t phi2)
    {
        return w;
    }
    //=============================================================================
    void _convection_impl(const FVField &phi,
                          const FVField &grad,
                          const FVBC &bc,
                          const Coefficient &vel,
                          la::MatSet lhs, la::VecSet rhs,
                          Interp auto &&interp, bool implicit)
    {
        // Finite volume space
        const auto V = phi.space();

        auto work = [&](const mesh::Mesh &mesh,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int facet_idx)
        {
            // Cells adjacent to the facet
            const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
            const auto &[cell_idx1, cell_idx2] = adjacent_cells;

            // Facet area and normal vector
            const real_t area = V->facet_area(facet_idx);
            const auto normal = V->facet_normal(facet_idx);

            // Facet geometric interpolation factor
            const real_t g = V->facet_interp_factor(facet_idx);

            // Compute facet normal flux
            real_t flux = 0.0;
            for (int i = 0; i < mesh.pdim(); i++)
            {
                const real_t facet_vel = g * vel(cell_idx1, i) + (1 - g) * vel(cell_idx2, i);
                flux += facet_vel * normal(i) * area;
            }

            // Cell values
            const real_t phi1 = phi(cell_idx1);
            const real_t phi2 = phi(cell_idx2);

            // Compute interpolation factor
            const real_t w = interp(flux, g, phi1, phi2);

            // Add contribution
            const std::array<int, 2> idxs = {cell_idx1, cell_idx2};
            if (implicit)
            {
                const std::array<real_t, 4> lhs_values = {w * flux, (1 - w) * flux,
                                                          -w * flux, -(1 - w) * flux};
                lhs(idxs, idxs, lhs_values);
            }
            else
            {
                const real_t phi_facet = w * phi1 + (1 - w) * phi2;
                const std::array<real_t, 2> rhs_values = {flux * phi_facet,
                                                          -flux * phi_facet};
                rhs(idxs, rhs_values);
            }
        };
        mesh::utils::for_all_facets(*V->mesh(), work, true, true);
    }
    //=============================================================================
    void convection(const FVField &phi,
                    const FVField &grad,
                    const FVBC &bc,
                    const Coefficient &vel,
                    la::MatSet lhs, la::VecSet rhs,
                    DifferencingScheme scheme, bool implicit)
    {
        SFEM_CHECK_SIZES(phi.space()->mesh()->pdim(), vel.n_comp());
        SFEM_CHECK_SIZES(phi.space()->mesh()->pdim(), grad.n_comp());

        switch (scheme)
        {
        case DifferencingScheme::upwind:
            _convection_impl(phi, grad, bc, vel, lhs, rhs, upwind, implicit);
            break;
        case DifferencingScheme::linear:
            _convection_impl(phi, grad, bc, vel, lhs, rhs, linear, implicit);
            break;
        default:
            SFEM_ERROR(std::format("Invalid differencing scheme: {}\n",
                                   static_cast<int>(scheme)));
            break;
        }
    }
}