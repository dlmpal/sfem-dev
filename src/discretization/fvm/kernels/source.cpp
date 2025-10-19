#include "source.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void add_source_term(const fvm::FVField &phi, la::VecSet vecset)
    {
        // Finite volume space
        const auto V = phi.space();

        // Store cell values
        std::vector<real_t> values(phi.n_comp());

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            std::array<int, 1> idx = {cell_idx};
            for (int i = 0; i < phi.n_comp(); i++)
            {
                values[i] = phi(cell_idx, i) * V->cell_volume(cell_idx);
            }
            vecset(idx, values);
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }
}