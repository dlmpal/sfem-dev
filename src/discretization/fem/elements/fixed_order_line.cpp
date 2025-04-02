#include "fixed_order.hpp"

namespace sfem::fem::fixed_order
{
    //=============================================================================
    void Line2::eval_shape(const std::array<real_t, 3> &pt, la::DenseMatrix &N) const
    {
        N(0, 0) = 0.5 * (1.0 - pt[0]);
        N(1, 0) = 0.5 * (1.0 + pt[0]);
    }
    //=============================================================================
    void Line2::eval_shape_grad(const std::array<real_t, 3> &pt, la::DenseMatrix &dNdxi) const
    {
        (void)pt;

        dNdxi(0, 0) = -0.5;
        dNdxi(1, 0) = 0.5;
    }

    //=============================================================================
    void Line3::eval_shape(const std::array<real_t, 3> &pt, la::DenseMatrix &N) const
    {
        N(0, 0) = -0.5 * pt[0] * (1.0 - pt[0]);
        N(1, 0) = 0.5 * pt[0] * (1 + pt[0]);
        N(2, 0) = 1.0 - pt[0] * pt[0];
    }
    //=============================================================================
    void Line3::eval_shape_grad(const std::array<real_t, 3> &pt, la::DenseMatrix &dNdxi) const
    {
        dNdxi(0, 0) = -0.5 + pt[0];
        dNdxi(1, 0) = 0.5 + pt[0];
        dNdxi(2, 0) = -2.0 * pt[0];
    }
    //=============================================================================
    void Line4::eval_shape(const std::array<real_t, 3> &pt, la::DenseMatrix &N) const
    {
        N(0, 0) = (0.333333333333333 - pt[0]) * (0.5625 * pt[0] - 0.5625) * (pt[0] + 0.333333333333333);
        N(1, 0) = (0.333333333333333 - pt[0]) * (-0.5625 * pt[0] - 0.5625) * (pt[0] + 0.333333333333333);
        N(2, 0) = (0.333333333333333 - pt[0]) * (1.6875 - 1.6875 * pt[0]) * (pt[0] + 1);
        N(3, 0) = (1.6875 - 1.6875 * pt[0]) * (pt[0] + 0.333333333333333) * (pt[0] + 1);
    }
    //=============================================================================
    void Line4::eval_shape_grad(const std::array<real_t, 3> &pt, la::DenseMatrix &dNdxi) const
    {
        dNdxi(0, 0) = -1.6875 * pt[0] * pt[0] + 1.125 * pt[0] + 0.0625;
        dNdxi(1, 0) = 1.6875 * pt[0] * pt[0] + 1.125 * pt[0] - 0.0625;
        dNdxi(2, 0) = 5.0625 * pt[0] * pt[0] - 1.125 * pt[0] - 1.6875;
        dNdxi(3, 0) = -5.0625 * pt[0] * pt[0] - 1.125 * pt[0] + 1.6875;
    }
}