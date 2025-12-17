#pragma once

#include "../elements/fe.hpp"
#include "../fe_field.hpp"

namespace sfem::fem::kernels
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