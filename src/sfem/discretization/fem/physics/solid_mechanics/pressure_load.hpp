#pragma once

#include <sfem/discretization/fem/core/elements/fe.hpp>
#include <sfem/discretization/fem/core/fe_field.hpp>

namespace sfem::fem::solid_mechanics
{
    class PressureLoad
    {
    public:
        PressureLoad(FEField U, Field &P, const mesh::Region &region);

        Field &P();
        const Field &P() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FEField U_;
        Field &P_;
        mesh::Region region_;
    };
}