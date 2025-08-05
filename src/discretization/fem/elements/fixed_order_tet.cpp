#include "fixed_order.hpp"

namespace sfem::fem::fixed_order
{
    //=============================================================================
    Tet4::Tet4()
        : NodalFiniteElement(mesh::CellType::tetrahedron, 1,
                             std::make_unique<quadrature::Tetrahedron>(4))
    {
    }
    //=============================================================================
    void Tet4::eval_shape(const std::array<real_t, 3> &pt, la::DenseMatrix &N) const
    {
        N(0, 0) = 1 - pt[0] - pt[1] - pt[2];
        N(1, 0) = pt[0];
        N(2, 0) = pt[1];
        N(3, 0) = pt[2];
    }
    //=============================================================================
    void Tet4::eval_shape_grad(const std::array<real_t, 3> &pt, la::DenseMatrix &dNdxi) const
    {
        (void)pt;

        dNdxi(0, 0) = -1.0;
        dNdxi(0, 1) = -1.0;
        dNdxi(0, 2) = -1.0;

        dNdxi(1, 0) = 1.0;
        dNdxi(1, 1) = 0.0;
        dNdxi(1, 2) = 0.0;

        dNdxi(2, 0) = 0.0;
        dNdxi(2, 1) = 1.0;
        dNdxi(2, 2) = 0.0;

        dNdxi(3, 0) = 0.0;
        dNdxi(3, 1) = 0.0;
        dNdxi(3, 2) = 1.0;
    }
    //=============================================================================
    Tet10::Tet10()
        : NodalFiniteElement(mesh::CellType::tetrahedron, 2,
                             std::make_unique<quadrature::Tetrahedron>(5))
    {
    }
    //=============================================================================
    void Tet10::eval_shape(const std::array<real_t, 3> &pt, la::DenseMatrix &N) const
    {
        real_t L1 = 1 - pt[0] - pt[1] - pt[2];
        real_t L2 = pt[0];
        real_t L3 = pt[1];
        real_t L4 = pt[2];

        // Corner nodes
        N(0, 0) = L1 * (2 * L1 - 1);
        N(1, 0) = L2 * (2 * L2 - 1);
        N(2, 0) = L3 * (2 * L3 - 1);
        N(3, 0) = L4 * (2 * L4 - 1);

        // Edge nodes
        N(4, 0) = 4 * L2 * L1;
        N(5, 0) = 4 * L2 * L3;
        N(6, 0) = 4 * L3 * L1;
        N(7, 0) = 4 * L4 * L1;
        N(8, 0) = 4 * L3 * L4;
        N(9, 0) = 4 * L2 * L4;
    }
    //=============================================================================
    void Tet10::eval_shape_grad(const std::array<real_t, 3> &pt, la::DenseMatrix &dNdxi) const
    {
        // Corner nodes
        dNdxi(0, 0) = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
        dNdxi(0, 1) = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
        dNdxi(0, 2) = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;

        dNdxi(1, 0) = 4 * pt[0] - 1;
        dNdxi(1, 1) = 0;
        dNdxi(1, 2) = 0;

        dNdxi(2, 0) = 0;
        dNdxi(2, 1) = 4 * pt[1] - 1;
        dNdxi(2, 2) = 0;

        dNdxi(3, 0) = 0;
        dNdxi(3, 1) = 0;
        dNdxi(3, 2) = 4 * pt[2] - 1;

        // Edge nodes
        dNdxi(4, 0) = -4 * (2 * pt[0] + pt[1] + pt[2] - 1);
        dNdxi(4, 1) = -4 * pt[0];
        dNdxi(4, 2) = -4 * pt[0];

        dNdxi(5, 0) = 4 * pt[1];
        dNdxi(5, 1) = 4 * pt[0];
        dNdxi(5, 2) = 0.0;

        dNdxi(6, 0) = -4 * pt[1];
        dNdxi(6, 1) = -4 * (pt[0] + 2 * pt[1] + pt[2] - 1);
        dNdxi(6, 2) = -4 * pt[1];

        dNdxi(7, 0) = -4 * pt[2];
        dNdxi(7, 1) = -4 * pt[2];
        dNdxi(7, 2) = -4 * (pt[0] + pt[1] + 2 * pt[2] - 1);

        dNdxi(8, 0) = 0;
        dNdxi(8, 1) = 4 * pt[2];
        dNdxi(8, 2) = 4 * pt[1];

        dNdxi(9, 0) = 4 * pt[2];
        dNdxi(9, 1) = 0;
        dNdxi(9, 2) = 4 * pt[0];
    }
}