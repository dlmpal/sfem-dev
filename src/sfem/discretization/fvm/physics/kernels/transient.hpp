#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/la/native/setval_utils.hpp>

namespace sfem::fvm
{
    class ImplicitEuler
    {
    public:
        ImplicitEuler(FVField phi, Field &C, real_t &dt);

        FVField &field();
        const FVField &field() const;

        Field &coeff();
        const Field &coeff() const;

        real_t &dt();
        real_t dt() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FVField phi_;

        Field &C_;

        real_t &dt_;
    };
}