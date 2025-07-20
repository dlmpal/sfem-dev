#include "diffusion.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void isotropic_diffusion(const FVFunction &phi,
                             const Coefficient &coeff,
                             const std::string &region,
                             la::MatSet lhs, la::VecSet rhs)
    {
        // (gradΦ)_f * S_f = (gradΦ)_f * (E_f + T_f)
        // (gradΦ)_f * S_f = E_f (Φ_F - Φ_C) / d_CF + gradΦ_f * T_f
        //  Distance between cell centers:
        //  e = (r_F - r_C) / ||r_F - r_C|| = d_CF / ||d_CF||

        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();

        // Integrate interior facets
        for (const auto &[facet, facet_idx] : mesh->region_facets(region))
        {
            // Integrate locally owned facets only
            if (mesh->topology()->entity_index_map(mesh->pdim() - 1)->is_ghost(facet_idx))
            {
                continue;
            }

            // Skip boundary facets
            if (V->is_boundary(facet_idx))
            {
                continue;
            }

            // Facet adjacent cells
            auto adjacent_cells = V->facet_adjacent_cells(facet_idx);

            // Facet area and normal vector
            const real_t area = V->facet_area(facet_idx);
            const auto normal = V->facet_normal(facet_idx);

            // Intercell distance and cell-facet normald istances
            const auto intercell_distance = V->intercell_distance(facet_idx);
            const real_t normal_distance = fabs(geo::inner(normal, intercell_distance));

            // Facet diffusion coefficient
            const real_t coeff_ = coeff(facet_idx, 0);

            // LHS matrix entry (orthogonal contribution)
            const real_t entry = coeff_ * area / normal_distance;
            std::array<real_t, 4> values = {entry, -entry, -entry, entry};
            lhs(adjacent_cells, adjacent_cells, values);

            // RHS vector entry (non-orthogonal contribution)
            /// @todo
        }
    }
    //=============================================================================
    void dirichlet_bc(const FVFunction &phi,
                      const Coefficient &coeff,
                      const std::string &region,
                      real_t value,
                      la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();

        for (auto [facet, facet_idx] : mesh->region_facets(region))
        {
            // Integrate locally owned facets only
            if (mesh->topology()->entity_index_map(mesh->pdim() - 1)->is_ghost(facet_idx))
            {
                continue;
            }

            const real_t area = V->facet_area(facet_idx);
            const real_t normal_distance = V->intercell_distance(facet_idx).mag();
            const real_t coeff = 2 * area / normal_distance;

            std::array<int, 1> idx = {V->facet_adjacent_cells(facet_idx).front()};
            std::array<real_t, 1> mat_value = {coeff};
            std::array<real_t, 1> vec_value = {coeff * value};

            lhs(idx, idx, mat_value);
            rhs(idx, vec_value);
        }
    }
}