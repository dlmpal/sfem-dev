#pragma once

#include "../elements/fe.hpp"
#include "../../function.hpp"

namespace sfem::fem::kernels
{
    class PressureLoad2D
    {
    public:
        PressureLoad2D(std::shared_ptr<Function> thick,
                       std::shared_ptr<Function> pressure);

        la::DenseMatrix operator()(int facet_idx,
                                   const fem::FEData &data,
                                   const geo::Vec3 &normal);

    private:
        std::shared_ptr<Function> thick_;
        std::shared_ptr<Function> pressure_;
    };

    // class PressureLoad3D
    // {
    // public:
    //     PressureLoad3D(real_t thick, real_t pressure_value);

    //     la::DenseMatrix operator()(int facet_idx,
    //                                const fem::FEData &data,
    //                                const geo::Vec3 &normal);

    // private:
    //     real_t pressure_value_;
    // };
}