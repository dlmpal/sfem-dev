#include "fv_function.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FVFunction::FVFunction(std::shared_ptr<const FVSpace> fv_space,
                           const std::vector<std::string> &components)
        : Function(fv_space->index_map(), components),
          fv_space_(fv_space)
    {
    }
    //=============================================================================
    std::shared_ptr<const FVSpace> FVFunction::space() const
    {
        return fv_space_;
    }
    //=============================================================================
    void gradient(const FVFunction &phi,
                  const std::string &region,
                  FVFunction &grad)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();
        const int n_comp = phi.n_comp();

        if (phi.n_comp() * dim != grad.n_comp())
        {
            SFEM_ERROR("");
        }

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

            // Facet's adjacent cells
            auto adjacent_cells = V->facet_adjacent_cells(facet_idx);

            // Facet area and normal vector
            const real_t area = V->facet_area(facet_idx);
            const auto normal = V->facet_normal(facet_idx);

            // Int{gradT}dV = Int{T}dS = Sum(Tf * Af * nf)

            for (int i = 0; i < n_comp; i++)
            {
                // Adjacent cell values
                const real_t phi_left = phi(adjacent_cells[0], i);
                const real_t phi_right = phi(adjacent_cells[1], i);

                // Weighted mean at facet
                const real_t phi_facet = 1;

                for (int j = 0; j < dim; j++)
                {
                    grad(adjacent_cells[0], i * n_comp + j) += phi_facet * area * normal(j);
                    grad(adjacent_cells[0], i * n_comp + j) -= phi_facet * area * normal(j);
                }
            }
        }
    }
}