#include "pressure_load.hpp"

namespace sfem::fem::kernels
{
    //=============================================================================
    PressureLoad2D::PressureLoad2D(std::shared_ptr<Function> thick,
                                   std::shared_ptr<Function> pressure)
        : thick_(thick), pressure_(pressure)
    {
    }
    //=============================================================================
    la::DenseMatrix PressureLoad2D::operator()(int facet_idx,
                                               const fem::FEData &data,
                                               const geo::Vec3 &normal)
    {
        // Facet coefficients
        const real_t thick = (*thick_)(facet_idx, 0);
        const real_t pressure_value = (*pressure_)(facet_idx, 0);

        la::DenseMatrix F(data.N.n_rows() * 2, 1, 0.0);
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            F(i * 2 + 0, 0) += -pressure_value * thick * data.N(i, 0) * normal.x();
            F(i * 2 + 1, 0) += -pressure_value * thick * data.N(i, 0) * normal.y();
        }
        return F;
    }
}