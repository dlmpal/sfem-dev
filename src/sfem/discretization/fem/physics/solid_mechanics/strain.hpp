#pragma once

#include <sfem/discretization/fem/core/elements/fe.hpp>
#include <sfem/discretization/fem/core/fe_field.hpp>

namespace sfem::fem::solid_mechanics
{
    class Strain
    {
    public:
        Strain(FEField U);

        FEField &U();
        const FEField &U() const;

        int n_strain() const;

        virtual void B_geo(const FEData &data, la::DenseMatrix &B) const = 0;
        virtual void B_mat(const FEData &data, la::DenseMatrix &B) const = 0;
        virtual void operator()(const FEData &data, la::DenseMatrix &e) const = 0;

    protected:
        FEField U_;
    };

    class SmallStrain : public Strain
    {
    public:
        SmallStrain(FEField U);
        void B_geo(const FEData &data, la::DenseMatrix &B) const override;
        void B_mat(const FEData &data, la::DenseMatrix &B) const override;
        void operator()(const FEData &data, la::DenseMatrix &e) const override;
    };
}