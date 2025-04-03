#include "fixed_order.hpp"

namespace sfem::fem::fixed_order
{
    //=============================================================================
    Tri3::Tri3()
        : NodalFiniteElement(mesh::CellType::triangle, 1,
                             std::make_unique<quadrature::Triangle>(3))
    {
    }
    //=============================================================================
    void Tri3::eval_shape(const std::array<real_t, 3> &pt, la::DenseMatrix &N) const
    {
        N(0, 0) = 1 - pt[0] - pt[1];
        N(1, 0) = pt[0];
        N(2, 0) = pt[1];
    }
    //=============================================================================
    void Tri3::eval_shape_grad(const std::array<real_t, 3> &pt, la::DenseMatrix &dNdxi) const
    {
        (void)pt;

        dNdxi(0, 0) = -1.0;
        dNdxi(0, 1) = -1.0;

        dNdxi(1, 0) = 1.0;
        dNdxi(1, 1) = 0.0;

        dNdxi(2, 0) = 0.0;
        dNdxi(2, 1) = 1.0;
    }
    //=============================================================================
    Tri6::Tri6()
        : NodalFiniteElement(mesh::CellType::triangle, 2,
                             std::make_unique<quadrature::Triangle>(4))
    {
    }
    //=============================================================================
    void Tri6::eval_shape(const std::array<real_t, 3> &pt, la::DenseMatrix &N) const
    {
        // Corner nodes
        N(0, 0) = (1 - pt[0] - pt[1]) * (1 - 2 * pt[0] - 2 * pt[1]);
        N(1, 0) = pt[0] * (2 * pt[0] - 1);
        N(2, 0) = pt[1] * (2 * pt[1] - 1);

        // Edge nodes
        N(3, 0) = 4 * pt[0] * (1 - pt[0] - pt[1]);
        N(4, 0) = 4 * pt[0] * pt[1];
        N(5, 0) = 4 * pt[1] * (1 - pt[0] - pt[1]);
    }
    //=============================================================================
    void Tri6::eval_shape_grad(const std::array<real_t, 3> &pt, la::DenseMatrix &dNdxi) const
    {
        // Corner nodes
        dNdxi(0, 0) = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
        dNdxi(0, 1) = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;

        dNdxi(1, 0) = 4.0 * pt[0] - 1.0;
        dNdxi(1, 1) = 0.0;

        dNdxi(2, 0) = 0.0;
        dNdxi(2, 1) = 4.0 * pt[1] - 1.0;

        // Edge nodes
        dNdxi(3, 0) = 4.0 - 8.0 * pt[0] - 4.0 * pt[1];
        dNdxi(3, 1) = -4.0 * pt[0];

        dNdxi(4, 0) = 4.0 * pt[1];
        dNdxi(4, 1) = 4.0 * pt[0];

        dNdxi(5, 0) = -4.0 * pt[1];
        dNdxi(5, 1) = 4.0 - 4.0 * pt[0] - 8.0 * pt[1];
    }
    //=============================================================================
    Tri10::Tri10()
        : NodalFiniteElement(mesh::CellType::triangle, 3,
                             std::make_unique<quadrature::Triangle>(6))
    {
    }
    //=============================================================================
    void Tri10::eval_shape(const std::array<real_t, 3> &pt, la::DenseMatrix &N) const
    {
        real_t L1 = 1 - pt[0] - pt[1];
        real_t L2 = pt[0];
        real_t L3 = pt[1];

        // Corner nodes
        N(0, 0) = 0.5 * (3 * L1 - 1) * (3 * L1 - 2) * L1;
        N(1, 0) = 0.5 * (3 * L2 - 1) * (3 * L2 - 2) * L2;
        N(2, 0) = 0.5 * (3 * L3 - 1) * (3 * L3 - 2) * L3;

        // Edge nodes
        N(3, 0) = 9 / 2 * L1 * L2 * (3 * L1 - 1);
        N(4, 0) = 9 / 2 * L1 * L2 * (3 * L2 - 1);

        N(5, 0) = 9 / 2 * L2 * L3 * (3 * L2 - 1);
        N(6, 0) = 9 / 2 * L2 * L3 * (3 * L3 - 1);

        N(7, 0) = 9 / 2 * L3 * L1 * (3 * L3 - 1);
        N(8, 0) = 9 / 2 * L3 * L1 * (3 * L1 - 1);

        // Internal node
        N(9, 0) = 27 * L1 * L2 * L3;
    }
    //=============================================================================
    void Tri10::eval_shape_grad(const std::array<real_t, 3> &pt, la::DenseMatrix &dNdxi) const
    {
        // Corner nodes
        dNdxi(0, 0) = -13.5 * pt[1] * pt[1] - 27.0 * pt[1] * pt[0] + 18.0 * pt[1] - 13.5 * pt[0] * pt[0] + 18.0 * pt[0] - 5.5;
        dNdxi(0, 1) = -13.5 * pt[1] * pt[1] - 27.0 * pt[1] * pt[0] + 18.0 * pt[1] - 13.5 * pt[0] * pt[0] + 18.0 * pt[0] - 5.5;

        dNdxi(1, 0) = 13.5 * pt[0] * pt[0] - 9.0 * pt[0] + 1.0;
        dNdxi(1, 1) = 0;

        dNdxi(2, 0) = 0;
        dNdxi(2, 1) = 13.5 * pt[1] * pt[1] - 9.0 * pt[1] + 1.0;

        // Edge nodes
        dNdxi(3, 0) = 13.5 * pt[1] * pt[1] + 54.0 * pt[1] * pt[0] - 22.5 * pt[1] + 40.5 * pt[0] * pt[0] - 45.0 * pt[0] + 9.0;
        dNdxi(3, 1) = pt[0] * (27.0 * pt[1] + 27.0 * pt[0] - 22.5);

        dNdxi(4, 0) = -27.0 * pt[1] * pt[0] + 4.5 * pt[1] - 40.5 * pt[0] * pt[0] + 36.0 * pt[0] - 4.5;
        dNdxi(4, 1) = pt[0] * (4.5 - 13.5 * pt[0]);

        dNdxi(5, 0) = pt[1] * (27.0 * pt[0] - 4.5);
        dNdxi(5, 1) = pt[0] * (13.5 * pt[0] - 4.5);

        dNdxi(6, 0) = pt[1] * (13.5 * pt[1] - 4.5);
        dNdxi(6, 1) = pt[0] * (27.0 * pt[1] - 4.5);

        dNdxi(7, 0) = pt[1] * (4.5 - 13.5 * pt[1]);
        dNdxi(7, 1) = -40.5 * pt[1] * pt[1] - 27.0 * pt[1] * pt[0] + 36.0 * pt[1] + 4.5 * pt[0] - 4.5;

        dNdxi(8, 0) = pt[1] * (27.0 * pt[1] + 27.0 * pt[0] - 22.5);
        dNdxi(8, 1) = 40.5 * pt[1] * pt[1] + 54.0 * pt[1] * pt[0] - 45.0 * pt[1] + 13.5 * pt[0] * pt[0] - 22.5 * pt[0] + 9.0;

        // Internal node
        dNdxi(9, 0) = 27 * pt[1] * (-pt[1] - 2 * pt[0] + 1);
        dNdxi(9, 1) = 27 * pt[0] * (-2 * pt[1] - pt[0] + 1);
    }
}