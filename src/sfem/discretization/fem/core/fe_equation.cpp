#include "fe_equation.hpp"
#include <sfem/discretization/fem/core/utils/la_utils.hpp>

namespace sfem::fem
{
    //=============================================================================
    Equation::Equation(FEField phi, std::shared_ptr<la::LinearSystem> Axb)
        : phi_(phi),
          bc_(phi_.space(), phi_.n_comp()),
          Axb_(Axb)
    {
        if (Axb_ == nullptr)
        {
            Axb_ = create_axb(phi_);
        }
    }
    //=============================================================================
    FEField &Equation::field()
    {
        return phi_;
    }
    //=============================================================================
    const FEField &Equation::field() const
    {
        return phi_;
    }
    //=============================================================================
    DirichletBC &Equation::bc()
    {
        return bc_;
    }
    //=============================================================================
    const DirichletBC &Equation::bc() const
    {
        return bc_;
    }
    //=============================================================================
    std::shared_ptr<la::LinearSystem> Equation::Axb()
    {
        return Axb_;
    }
    //=============================================================================
    std::shared_ptr<const la::LinearSystem> Equation::Axb() const
    {
        return Axb_;
    }
    //=============================================================================
    Equation &Equation::add_kernel(const FEKernel &kernel)
    {
        kernels_.push_back(kernel);
        return *this;
    }
    //=============================================================================
    void Equation::clear_kernels()
    {
        kernels_.clear();
    }
    //=============================================================================
    void Equation::assemble()
    {
        Axb_->reset();

        for (const FEKernel &kernel : kernels_)
        {
            kernel(Axb_->lhs(), Axb_->rhs());
        }

        Axb_->assemble();
    }
    //=============================================================================
    void Equation::apply_dirichlet_bc()
    {
        const auto [cdof, cvalues] = bc_.get_dofs_values();
        Axb_->eliminate_dofs(cdof, cvalues);

        auto &values = phi_.dof_values().values();
        for (std::size_t i = 0; i < cdof.size(); i++)
        {
            values[cdof[i]] = cvalues[i];
        }
    }
    //=============================================================================
    void Equation::solve()
    {
        Axb_->solve(phi_.dof_values());
        phi_.dof_values().update_ghosts();
    }
}