#pragma once

#include <sfem/discretization/fem/core/elements/fe.hpp>
#include <sfem/discretization/fem/core/fe_field.hpp>

namespace sfem::fem
{
    class Diffusion
    {
    public:
        Diffusion(FEField phi, Field &D);

        FEField field() const;

        Field &D();
        const Field &D() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FEField phi_;
        Field &D_;
    };
}