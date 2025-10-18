#include "diffusion.hpp"
#include "laplacian.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void diffusion(const FVField &phi, const FVField &grad,
                   const FVBC &bc, real_t dt,
                   const Coefficient &cell_coeff,
                   const Coefficient &facet_coeff,
                   la::MatSet lhs, la::VecSet rhs)
    {
        // Assemble Laplacian
        laplacian(phi, grad, bc, facet_coeff, lhs, rhs);

        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();

        const real_t dt_inv = 1.0 / dt;
        for (const auto &region : mesh->regions())
        {
            if (region.dim() < dim)
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

                const real_t coeff = cell_coeff(cell_idx, 0);
                const real_t vol = V->cell_volume(cell_idx);
                std::array<int, 1> idx = {cell_idx};
                std::array<real_t, 1> lhs_value = {coeff * vol * dt_inv};
                std::array<real_t, 1> rhs_value = {coeff * vol * dt_inv * phi(cell_idx)};
                lhs(idx, idx, lhs_value);
                rhs(idx, rhs_value);
            }
        }
    }
}