#include "pressure_load.hpp"

namespace sfem::fem::kernels
{
    //=============================================================================
    PressureLoad2D::PressureLoad2D(real_t thick, real_t pressure_value)
        : thick_(thick), pressure_value_(pressure_value)
    {
    }
    //=============================================================================
    la::DenseMatrix PressureLoad2D::operator()(int facet_idx,
                                               const fem::FEData &data,
                                               const geo::Vec3 &normal)
    {
        la::DenseMatrix F(data.N.n_rows() * 2, 1, 0.0);
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            F(i * 2 + 0, 0) += -pressure_value_ * thick_ * data.N(i, 0) * normal.x();
            F(i * 2 + 1, 0) += -pressure_value_ * thick_ * data.N(i, 0) * normal.y();
        }
        return F;
    }
}