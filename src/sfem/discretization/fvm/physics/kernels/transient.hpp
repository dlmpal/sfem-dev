#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/la/native/setval_utils.hpp>

namespace sfem::fvm
{
    class ImplicitEuler
    {
    public:
        ImplicitEuler(FVField phi, IField &C, real_t &dt);

        FVField &field();
        const FVField &field() const;

        IField &coeff();
        const IField &coeff() const;

        real_t &dt();
        real_t dt() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FVField phi_;

        IField &C_;

        real_t &dt_;
    };
}