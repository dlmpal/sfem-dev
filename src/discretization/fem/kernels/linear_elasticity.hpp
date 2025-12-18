#pragma once

#include "../elements/fe.hpp"
#include "../fe_field.hpp"

namespace sfem::fem::kernels::elasticity
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

    /// @brief Todo make abstract
    class Strain
    {
    public:
        Strain(FEField U);
        int n_strain(int dim) const;
        void B_geo(const FEData &data, la::DenseMatrix &B) const;
        void B_mat(const FEData &data, la::DenseMatrix &B) const;
        void operator()(const FEData &data, la::DenseMatrix &e) const;

    protected:
        FEField U_;
    };

    class Stress
    {
    public:
        Stress(FEField U, Strain &strain, ElasticityConstitutive &constitutive);
        /// @todo add accessors
        void operator()(const FEData &data, la::DenseMatrix &e) const;

    private:
        FEField U_;
        Strain &strain_;
        ElasticityConstitutive &constitutive_;
    };

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