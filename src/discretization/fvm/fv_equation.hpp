#pragma once

#include "fv_field.hpp"
#include "../../la/native/linear_system.hpp"

namespace sfem::fvm
{
    using FVKernel = std::function<void(la::MatSet, la::VecSet)>;

    class Equation
    {
    public:
        Equation(FVField phi, std::shared_ptr<la::LinearSystem> Axb = nullptr);

        FVField field() const;

        std::shared_ptr<la::LinearSystem> Axb();

        std::shared_ptr<const la::LinearSystem> Axb() const;

        const la::Vector &diag() const;

        Equation &add_kernel(const FVKernel &kernel);

        void clear_kernels();

        void assemble();

        void apply_relaxation(real_t alpha);

        void solve();

    protected:
        FVField phi_;

        std::vector<FVKernel> kernels_;

        std::shared_ptr<la::LinearSystem> Axb_;

        la::Vector diag_;
    };
}