#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/la/native/setval_utils.hpp>

namespace sfem::fvm
{
    class Source
    {
    public:
        using SourceFunc = std::function<void(const FVField &, int, std::span<real_t>)>;

        Source(FVField phi, SourceFunc func);

        FVField &field();
        const FVField &field() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FVField phi_;

        SourceFunc func_;
    };
}