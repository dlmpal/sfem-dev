#include "transient.hpp"
#include "../../../la/native/setval_utils.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    ImplicitEuler::ImplicitEuler(FVField phi, IField &C, real_t &dt)
        : phi_(phi),
          C_(C),
          dt_(dt)
    {
    }
    //=============================================================================
    FVField ImplicitEuler::field() const
    {
        return phi_;
    }
    //=============================================================================
    IField &ImplicitEuler::coeff()
    {
        return C_;
    }
    //=============================================================================
    const IField &ImplicitEuler::coeff() const
    {
        return C_;
    }
    //=============================================================================
    real_t &ImplicitEuler::dt()
    {
        return dt_;
    }
    //=============================================================================
    real_t ImplicitEuler::dt() const
    {
        return dt_;
    }
    //=============================================================================
    void ImplicitEuler::operator()(la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = phi_.space();

        // Invert once
        const real_t dt_inv = 1.0 / dt_;

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            const real_t vol = V->cell_volume(cell_idx);
            const std::array<int, 1> idx = {cell_idx};
            const std::array<real_t, 1> lhs_value = {C_.cell_value(cell_idx) * vol * dt_inv};
            const std::array<real_t, 1> rhs_value = {lhs_value[0] * phi_.cell_value(cell_idx)};
            lhs(idx, idx, lhs_value);
            rhs(idx, rhs_value);
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }
}