#pragma once

#include "../elements/fe.hpp"
#include "../../function.hpp"

namespace sfem::fem::kernels
{
    class LinearElasticity2D
    {
    public:
        LinearElasticity2D(std::shared_ptr<Function> E,
                           std::shared_ptr<Function> nu,
                           std::shared_ptr<Function> thick);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data) const;

    private:
        std::shared_ptr<Function> E_;
        std::shared_ptr<Function> nu_;
        std::shared_ptr<Function> thick_;
    };

    class LinearElasticity3D
    {
    public:
        LinearElasticity3D(std::shared_ptr<Function> E,
                           std::shared_ptr<Function> nu);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data) const;

    private:
        std::shared_ptr<Function> E_;
        std::shared_ptr<Function> nu_;
    };
}