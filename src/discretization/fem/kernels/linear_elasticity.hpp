#pragma once

#include "../elements/fe.hpp"
#include "../fe_field.hpp"

namespace sfem::fem::kernels::elasticity
{
    class LinearElasticIsotropic
    {
    public:
        LinearElasticIsotropic(Field &E, Field &nu, Field &rho);
        ~LinearElasticIsotropic() = default;

        Field &E();
        const Field &E() const;

        Field &nu();
        const Field &nu() const;

        Field &rho();
        const Field &rho() const;

        virtual int dim() const = 0;
        virtual int n_strain() const = 0;

        virtual void tangent(int cell_idx,
                             mesh::CellType cell_type,
                             const std::array<real_t, 3> &pt,
                             la::DenseMatrix &D) const = 0;

        virtual void stress(int cell_idx,
                            mesh::CellType cell_type,
                            const std::array<real_t, 3> &pt,
                            const la::DenseMatrix &strain,
                            la::DenseMatrix &stress) const = 0;

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

        void tangent(int cell_idx,
                     mesh::CellType cell_type,
                     const std::array<real_t, 3> &pt,
                     la::DenseMatrix &D) const override;

        void stress(int cell_idx,
                    mesh::CellType cell_type,
                    const std::array<real_t, 3> &pt,
                    const la::DenseMatrix &strain,
                    la::DenseMatrix &stress) const override;
    };

    class LinearElastic3D : public LinearElasticIsotropic
    {
    public:
        LinearElastic3D(Field &E, Field &nu, Field &rho);

        int dim() const override;
        int n_strain() const override;

        void tangent(int cell_idx,
                     mesh::CellType cell_type,
                     const std::array<real_t, 3> &pt,
                     la::DenseMatrix &D) const override;

        void stress(int cell_idx,
                    mesh::CellType cell_type,
                    const std::array<real_t, 3> &pt,
                    const la::DenseMatrix &strain,
                    la::DenseMatrix &stress) const override;
    };

    class LinearElasticity
    {
    public:
        LinearElasticity(FEField U, LinearElasticIsotropic &constitutive,
                         const std::array<real_t, 3> &g = {});

        void operator()(la::MatSet lhs, la::VecSet rhs);

        void strain_displacement(const FEData &data, la::DenseMatrix &B) const;

    private:
        FEField U_;
        LinearElasticIsotropic &constitutive_;
        std::array<real_t, 3> g_;
    };

    class PressureLoad
    {
    public:
        PressureLoad(FEField U, ConstantField &P, const mesh::Region &region);

        ConstantField &P();
        const ConstantField &P() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FEField U_;
        ConstantField &P_;
        mesh::Region region_;
    };
}