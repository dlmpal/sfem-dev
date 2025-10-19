#include "diffusion.hpp"
#include "laplacian.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void diffusion(const FVField &phi, const FVField &grad,
                   const FVBC &bc, real_t dt,
                   const Coefficient &cell_coeff,
                   const Coefficient &facet_coeff,
                   la::MatSet lhs, la::VecSet rhs)
    {
        // Add Laplacian contribution
        laplacian(phi, grad, bc, facet_coeff, lhs, rhs);

        // Finite volume space
        const auto V = phi.space();

        // Add transient term contribution
        const real_t dt_inv = 1.0 / dt;
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            const real_t coeff = cell_coeff(cell_idx, 0);
            const real_t vol = V->cell_volume(cell_idx);
            std::array<int, 1> idx = {cell_idx};
            std::array<real_t, 1> lhs_value = {coeff * vol * dt_inv};
            std::array<real_t, 1> rhs_value = {coeff * vol * dt_inv * phi(cell_idx)};
            lhs(idx, idx, lhs_value);
            rhs(idx, rhs_value);
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }
}