#pragma once

#include "../elements/fe.hpp"
#include "../../coefficient.hpp"

namespace sfem::fem::kernels
{
    class MassND
    {
    public:
        MassND(int n_comp, std::shared_ptr<const Coefficient> coeff);

        la::DenseMatrix operator()(int cell_idx, const FEData &data);

    private:
        int n_comp_;
        std::shared_ptr<const Coefficient> coeff_;
    };
}