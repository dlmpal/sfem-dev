#include "laplacian.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void laplacian(const FVField &phi,
                   const FVField &grad,
                   const FVBC &bc,
                   const Coefficient &coeff,
                   la::MatSet lhs, la::VecSet rhs)
    {
        if (phi.n_comp() > 1)
        {
            SFEM_ERROR(std::format("Only scalar functions are supported, (n_comp={} > 1)\n", phi.n_comp()));
        }

        // Finite volume space
        const auto V = phi.space();

        auto work = [&](const mesh::Mesh &mesh,
                        const mesh::Region &region,
                        const mesh::Cell &,
                        int facet_idx)
        {
            // Cells adjacent to the facet
            const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
            const auto &[cell_idx1, cell_idx2] = adjacent_cells;

            // Facet area and normal vector
            const real_t area = V->facet_area(facet_idx);
            const auto normal = V->facet_normal(facet_idx);

            // Distance of adjacent cell midpoints
            const auto d12 = V->intercell_distance(facet_idx);

            // Boundary facets
            if (cell_idx1 == cell_idx2)
            {
                if (bc.types_.contains(region.name()))
                {
                    const std::array<int, 1> idx = {cell_idx1};
                    std::array<real_t, 1> lhs_value{};
                    std::array<real_t, 1> rhs_value{};

                    if (bc.types_.at(region.name()) == BCType::dirichlet)
                    {
                        lhs_value[0] = 2.0 * coeff(cell_idx1, 0) * area / d12.mag();
                        rhs_value[0] = lhs_value[0] * bc.values_.at(facet_idx)[0];
                    }
                    else if (bc.types_.at(region.name()) == BCType::neumann)
                    {
                        lhs_value[0] = 0.0;
                        rhs_value[0] = coeff(cell_idx1, 0) * area * bc.values_.at(facet_idx)[0];
                    }
                    else if (bc.types_.at(region.name()) == BCType::robin)
                    {
                        const real_t phi_inf = bc.values_.at(facet_idx)[0];
                        const real_t h_inf = bc.values_.at(facet_idx)[1];
                        const real_t d = d12.mag();
                        lhs_value[0] = (h_inf * coeff(cell_idx1, 0) / d) /
                                       (h_inf + coeff(cell_idx1, 0) / d) * area;
                        rhs_value[0] = lhs_value[0] * phi_inf;
                    }

                    lhs(idx, idx, lhs_value);
                    rhs(idx, rhs_value);
                }
            }
            // Internal facets
            else
            {
                // Facet geometric interpolation factor
                const real_t g = V->facet_interp_factor(facet_idx);

                // Compute facet value from adjacent cell values
                const real_t coeff1 = coeff(cell_idx1, 0);
                const real_t coeff2 = coeff(cell_idx2, 0);
                const real_t coeff_facet = g * coeff1 + (1 - g) * coeff2;

                // Facet normal decomposition vectors
                const geo::Vec3 delta = (d12 / geo::inner(normal, d12));
                const geo::Vec3 kappa = (normal - delta);

                // LHS - orthogonal contribution
                const real_t value = coeff_facet * area / d12.mag();
                std::array<real_t, 4> lhs_values = {value, -value, -value, value};
                lhs(adjacent_cells, adjacent_cells, lhs_values);

                // RHS - non-orthogonal correction
                std::array<real_t, 2> rhs_values{};
                for (int i = 0; i < mesh.pdim(); i++)
                {
                    const real_t grad1 = grad(cell_idx1, i);
                    const real_t grad2 = grad(cell_idx2, i);
                    const real_t grad_facet = g * grad1 + (1 - g) * grad2;
                    rhs_values[0] += coeff_facet * grad_facet * kappa(i) * area;
                    rhs_values[1] -= coeff_facet * grad_facet * kappa(i) * area;
                }
                rhs(adjacent_cells, rhs_values);
            }
        };
        mesh::utils::for_all_facets(*V->mesh(), work);
    }
}