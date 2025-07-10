#pragma once

#include "../elements/fe.hpp"

namespace sfem::fem::kernels
{
    class Diffusion2D
    {
    public:
        Diffusion2D(real_t coeff);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data);

    private:
        real_t coeff_;
    };

    class Diffusion3D
    {
    public:
        Diffusion3D(real_t coeff);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data);

    private:
        real_t coeff_;
    };
}