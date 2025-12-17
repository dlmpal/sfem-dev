#pragma once

#include "../elements/fe.hpp"
#include "../fe_field.hpp"

namespace sfem::fem::kernels
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