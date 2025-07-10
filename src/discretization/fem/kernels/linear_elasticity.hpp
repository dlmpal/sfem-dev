#pragma once

#include "../elements/fe.hpp"

namespace sfem::fem::kernels
{
    class LinearElasticity2D
    {
    public:
        LinearElasticity2D(real_t E, real_t nu, real_t thick);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data) const;

    private:
        real_t E_;
        real_t nu_;
        real_t thick_;
    };

    // class LinearElasticity3D
    // {
    // public:
    //     LinearElasticity3D(real_t E, real_t nu);
    //     la::DenseMatrix operator()(int cell_idx, const fem::FEData &data) const;

    // private:
    //     real_t E_;
    //     real_t nu_;
    // };
}