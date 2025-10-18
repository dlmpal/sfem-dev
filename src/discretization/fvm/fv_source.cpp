#include "fv_source.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void add_source_term(const fvm::FVField &phi, la::VecSet vecset)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();

        // Loop over all cells
        std::vector<real_t> cell_values(phi.n_comp());
        for (const auto &region : mesh->regions())
        {
            // Skip boundary regions
            if (region.dim() < mesh->pdim())
            {
                continue;
            }

            for (const auto &[cell, cell_idx] : mesh->region_cells(region.name()))
            {
                // Skip ghost cells
                if (mesh->topology()->entity_index_map(dim)->is_ghost(cell_idx))
                {
                    continue;
                }

                std::array<int, 1> idx = {cell_idx};
                for (int i = 0; i < phi.n_comp(); i++)
                {
                    cell_values[i] = phi(cell_idx, i) * V->cell_volume(cell_idx);
                }
                vecset(idx, cell_values);
            }
        }
    }
}