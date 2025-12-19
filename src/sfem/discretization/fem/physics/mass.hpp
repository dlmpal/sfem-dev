#pragma once

#include <sfem/discretization/fem/core/elements/fe.hpp>
#include <sfem/discretization/fem/core/fe_field.hpp>

namespace sfem::fem
{
    class MassND
    {
    public:
        MassND(FEField phi, Field &C);

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FEField phi_;
        Field &C_;
    };
}