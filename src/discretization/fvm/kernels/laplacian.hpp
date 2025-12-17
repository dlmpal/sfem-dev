#pragma once

#include "../fv_field.hpp"
#include "../../../la/native/setval_utils.hpp"

namespace sfem::fvm
{
    class Laplacian
    {
    public:
        Laplacian(FVField phi, IField &D);

        FVField field() const;

        void operator()(la::MatSet lhs, la::VecSet rhs);

    private:
        FVField phi_;
        IField &D_;
    };
}