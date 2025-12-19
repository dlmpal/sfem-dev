#include "source.hpp"
#include <sfem/la/native/setval_utils.hpp>
#include <sfem/mesh/utils/loop_utils.hpp>

namespace sfem::fvm
{
    //=============================================================================
    Source::Source(FVField phi, SourceFunc func)
        : phi_(phi),
          func_(func)
    {
    }
    //=============================================================================
    FVField Source::field() const
    {
        return phi_;
    }
    //=============================================================================
    void Source::operator()(la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = phi_.space();

        // Store cell values
        std::vector<real_t> values(phi_.n_comp());

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            std::array<int, 1> idx = {cell_idx};
            func_(phi_, cell_idx, values);
            for (auto &value : values)
            {
                value *= V->cell_volume(cell_idx);
            }
            rhs(idx, values);
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }
}