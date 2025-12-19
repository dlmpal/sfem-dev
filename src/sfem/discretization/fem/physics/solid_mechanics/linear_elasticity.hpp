#pragma once

#include <sfem/discretization/fem/physics/solid_mechanics/constitutive.hpp>
#include <sfem/discretization/fem/physics/solid_mechanics/strain.hpp>

namespace sfem::fem::solid_mechanics
{
    class LinearElasticity
    {
    public:
        LinearElasticity(FEField U, Strain &strain,
                         LinearElasticIsotropic &constitutive,
                         const std::array<real_t, 3> &g = {});

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FEField U_;
        Strain &strain_;
        LinearElasticIsotropic &constitutive_;
        std::array<real_t, 3> g_;
    };
}