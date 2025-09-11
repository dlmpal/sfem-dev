#include "diffusion.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void diffusion(const FVFunction &phi,
                   const FVFunction &grad,
                   const FVBC &bc,
                   const Coefficient &coeff,
                   la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();

        // Integrate facets
        for (const auto &region : mesh->regions())
        {
            for (const auto &[facet, facet_idx] : mesh->region_facets(region.name()))
            {
                // Integrate locally owned facets only
                if (mesh->topology()->entity_index_map(mesh->pdim() - 1)->is_ghost(facet_idx))
                {
                    continue;
                }

                // Facet adjacent cells
                const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
                const auto &[cell_idx1, cell_idx2] = adjacent_cells;

                // Facet area and normal vector
                const real_t area = V->facet_area(facet_idx);
                const auto normal = V->facet_normal(facet_idx);

                // Intercell and normal distance
                const auto d12 = V->intercell_distance(facet_idx);
                const real_t d12_normal = geo::inner(normal, d12);

                //
                geo::Vec3 delta = d12 / d12_normal;

                //
                geo::Vec3 kappa = normal - delta;

                if (adjacent_cells[0] == adjacent_cells[1])
                {
                    if (bc.types_.contains(region.name()))
                    {
                        if (bc.types_.at(region.name()) == BCType::dirichlet)
                        {
                            const real_t coeff_facet = coeff(cell_idx1, 0);
                            const real_t entry = 2 * coeff_facet * area / d12.mag();

                            std::array<int, 1> idx = {cell_idx1};
                            std::array<real_t, 1> lhs_value = {entry};
                            std::array<real_t, 1> rhs_value = {entry * bc.values_.at(facet_idx)};

                            lhs(idx, idx, lhs_value);
                            rhs(idx, rhs_value);
                        }
                    }
                }
                else
                {
                    // Facet diffusion coefficient
                    const real_t coeff1 = coeff(cell_idx1, 0);
                    const real_t coeff2 = coeff(cell_idx2, 0);
                    const real_t coeff_facet = V->compute_facet_value(facet_idx,
                                                                      coeff1, coeff2);

                    // LHS - orthogonal contribution
                    const real_t value = coeff_facet * area / d12_normal;
                    std::array<real_t, 4> lhs_values = {value, -value, -value, value};
                    lhs(adjacent_cells, adjacent_cells, lhs_values);

                    // RHS - non-orthogonal correction
                    std::array<real_t, 2> rhs_values{};
                    for (int i = 0; i < dim; i++)
                    {
                        const real_t grad1 = grad(cell_idx1, i);
                        const real_t grad2 = grad(cell_idx2, i);
                        const real_t grad_facet = V->compute_facet_value(facet_idx,
                                                                         grad1, grad2);
                        rhs_values[0] += coeff_facet * grad_facet * kappa(i) * area;
                        rhs_values[1] -= coeff_facet * grad_facet * kappa(i) * area;
                    }
                    rhs(adjacent_cells, rhs_values);
                }
            }
        }
    }
}