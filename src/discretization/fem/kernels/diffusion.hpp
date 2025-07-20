#pragma once

#include "../elements/fe.hpp"
#include "../../coefficient.hpp"

namespace sfem::fem::kernels
{
    class Diffusion2D
    {
    public:
        Diffusion2D(std::shared_ptr<const Coefficient> coeff);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data);

    private:
        std::shared_ptr<const Coefficient> coeff_;
    };

    class Diffusion3D
    {
    public:
        Diffusion3D(std::shared_ptr<const Coefficient> coeff);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data);

    private:
        std::shared_ptr<const Coefficient> coeff_;
    };
}