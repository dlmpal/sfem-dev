#pragma once

#include <sfem/discretization/fem/core/fe_field.hpp>
#include <sfem/discretization/fem/core/dirichlet_bc.hpp>
#include <sfem/la/native/linear_system.hpp>

namespace sfem::fem
{
    using FEKernel = std::function<void(la::MatSet, la::VecSet)>;

    class Equation
    {
    public:
        Equation(FEField phi, std::shared_ptr<la::LinearSystem> Axb = nullptr);

        FEField &field();
        const FEField &field() const;

        DirichletBC &bc();
        const DirichletBC &bc() const;

        std::shared_ptr<la::LinearSystem> Axb();
        std::shared_ptr<const la::LinearSystem> Axb() const;

        Equation &add_kernel(const FEKernel &kernel);

        void clear_kernels();

        void assemble();

        void apply_dirichlet_bc();

        void solve();

    protected:
        FEField phi_;

        DirichletBC bc_;

        std::shared_ptr<la::LinearSystem> Axb_;

        std::vector<FEKernel> kernels_;
    };
}