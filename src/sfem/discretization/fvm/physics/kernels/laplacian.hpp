#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/la/native/setval_utils.hpp>

namespace sfem::fvm
{
    class Laplacian
    {
    public:
        Laplacian(FVField phi, IField &D);

        FVField &field();
        const FVField &field() const;

        IField &D();
        const IField &D() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FVField phi_;
        IField &D_;
    };
}