#include "fv_equation.hpp"
#include "utils/la_utils.hpp"
#include "../../la/native/setval_utils.hpp"
#include "../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    Equation::Equation(FVField phi, std::shared_ptr<la::LinearSystem> Axb)
        : phi_(phi),
          Axb_(Axb),
          diag_(create_vec(phi))
    {
        if (Axb_ == nullptr)
        {
            Axb_ = create_axb(phi_);
        }
    }
    //=============================================================================
    FVField Equation::field() const
    {
        return phi_;
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
    const la::Vector &Equation::diag() const
    {
        return diag_;
    }
    //=============================================================================
    Equation &Equation::add_kernel(const FVKernel &kernel)
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

        for (const FVKernel &kernel : kernels_)
        {
            kernel(Axb_->lhs(), Axb_->rhs());
        }

        Axb_->assemble();

        Axb_->diagonal(diag_);
    }
    //=============================================================================
    void Equation::apply_relaxation(real_t alpha)
    {
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            diag_(cell_idx) *= (1 - alpha) / alpha * phi_.cell_value(cell_idx);
        };
        mesh::utils::for_all_cells(*phi_.space()->mesh(), work);

        Axb_->rhs_axpy(1.0, diag_);

        Axb_->scale_diagonal(1.0 / alpha);

        Axb_->diagonal(diag_);
    }
    //=============================================================================
    void Equation::solve()
    {
        Axb_->solve(phi_.values());
        phi_.values().update_ghosts();
        phi_.update_gradient();
    }
}