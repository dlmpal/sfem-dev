#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/la/native/setval_utils.hpp>

namespace sfem::fvm
{
    class Laplacian
    {
    public:
        Laplacian(FVField phi, Field &D);

        FVField &field();
        const FVField &field() const;

        Field &D();
        const Field &D() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FVField phi_;
        Field &D_;
    };
}