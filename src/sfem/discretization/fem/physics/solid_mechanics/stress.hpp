#pragma once

#include <sfem/discretization/fem/physics/solid_mechanics/constitutive.hpp>
#include <sfem/discretization/fem/physics/solid_mechanics/strain.hpp>

namespace sfem::fem::solid_mechanics
{
    class Stress
    {
    public:
        Stress(FEField U, Strain &strain, ElasticityConstitutive &constitutive);

        FEField &U();
        const FEField &U() const;

        Strain &strain();
        const Strain &strain() const;

        ElasticityConstitutive &constitutive();
        const ElasticityConstitutive &constitutive() const;

        void operator()(const FEData &data, la::DenseMatrix &e) const;

    private:
        FEField U_;
        Strain &strain_;
        ElasticityConstitutive &constitutive_;
    };
}