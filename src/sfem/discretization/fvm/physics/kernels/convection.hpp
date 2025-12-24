#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/la/native/setval_utils.hpp>

namespace sfem::fvm
{
    class Convection
    {
    public:
        Convection(FVField phi, std::vector<real_t> &flux);

        FVField &field();
        const FVField &field() const;

        std::vector<real_t> &flux();
        const std::vector<real_t> &flux() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FVField phi_;

        /// @todo
        std::vector<real_t> &flux_;
    };
}