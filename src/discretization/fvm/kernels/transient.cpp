#include "transient.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void ddt(const FVField &phi, const Coefficient &coeff,
             real_t dt, la::MatSet lhs, la::VecSet rhs)
    {
        // Finite volume space
        const auto V = phi.space();

        const real_t dt_inv = 1.0 / dt;

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            const real_t vol = V->cell_volume(cell_idx);
            const std::array<int, 1> idx = {cell_idx};
            const std::array<real_t, 1> lhs_value = {coeff(cell_idx, 0) * vol * dt_inv};
            const std::array<real_t, 1> rhs_value = {lhs_value[0] * phi(cell_idx)};
            lhs(idx, idx, lhs_value);
            rhs(idx, rhs_value);
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }
}