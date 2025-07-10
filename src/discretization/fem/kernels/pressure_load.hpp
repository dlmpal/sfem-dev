#pragma once

#include "../elements/fe.hpp"

namespace sfem::fem::kernels
{
    class PressureLoad2D
    {
    public:
        PressureLoad2D(real_t thick, real_t pressure_value);

        la::DenseMatrix operator()(int facet_idx,
                                   const fem::FEData &data,
                                   const geo::Vec3 &normal);

    private:
        real_t thick_;
        real_t pressure_value_;
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