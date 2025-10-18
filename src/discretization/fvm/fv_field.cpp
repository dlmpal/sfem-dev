#include "fv_field.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FVField::FVField(std::shared_ptr<const FVSpace> fv_space,
                     const std::vector<std::string> &components)
        : Field(fv_space->index_map(), components),
          fv_space_(fv_space)
    {
    }
    //=============================================================================
    std::shared_ptr<const FVSpace> FVField::space() const
    {
        return fv_space_;
    }
    //=============================================================================
    void eval_field(FVField &phi, FieldFunction func, real_t time)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();

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
                if (mesh->topology()->entity_index_map(mesh->pdim())->is_ghost(cell_idx))
                {
                    continue;
                }

                //
                for (int i = 0; i < phi.n_comp(); i++)
                {
                    cell_values[i] = phi(cell_idx, i);
                }

                //
                func(V->cell_midpoint(cell_idx), cell_values, time);

                //
                for (int i = 0; i < phi.n_comp(); i++)
                {
                    phi(cell_idx, i) = cell_values[i];
                }
            }
        }

        //
        phi.update_ghosts();
    }
}