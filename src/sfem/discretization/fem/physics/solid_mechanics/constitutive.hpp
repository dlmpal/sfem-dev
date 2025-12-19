#pragma once

#include <sfem/discretization/fem/core/elements/fe.hpp>
#include <sfem/discretization/fem/core/fe_field.hpp>

namespace sfem::fem::solid_mechanics
{
    class ElasticityConstitutive
    {
    public:
        virtual ~ElasticityConstitutive() = default;

        virtual int dim() const = 0;
        virtual int n_strain() const = 0;

        virtual void stress(const FEData &data,
                            const la::DenseMatrix &strain,
                            la::DenseMatrix &stress) const = 0;

        virtual void tangent(const FEData &data,
                             const la::DenseMatrix &strain,
                             la::DenseMatrix &D) const = 0;
    };

    class LinearElasticIsotropic : public ElasticityConstitutive
    {
    public:
        LinearElasticIsotropic(Field &E, Field &nu, Field &rho);

        Field &E();
        const Field &E() const;

        Field &nu();
        const Field &nu() const;

        Field &rho();
        const Field &rho() const;

    protected:
        Field &E_;
        Field &nu_;
        Field &rho_;
    };

    class LinearElasticPlaneStress : public LinearElasticIsotropic
    {
    public:
        LinearElasticPlaneStress(Field &E, Field &nu, Field &rho);

        int dim() const override;
        int n_strain() const override;

        void stress(const FEData &data,
                    const la::DenseMatrix &strain,
                    la::DenseMatrix &stress) const override;

        void tangent(const FEData &data,
                     const la::DenseMatrix &strain,
                     la::DenseMatrix &D) const override;
    };

    class LinearElastic3D : public LinearElasticIsotropic
    {
    public:
        LinearElastic3D(Field &E, Field &nu, Field &rho);

        int dim() const override;
        int n_strain() const override;

        void stress(const FEData &data,
                    const la::DenseMatrix &strain,
                    la::DenseMatrix &stress) const override;

        void tangent(const FEData &data,
                     const la::DenseMatrix &strain,
                     la::DenseMatrix &D) const override;
    };
}