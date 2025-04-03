#include "fixed_order.hpp"

namespace sfem::fem::fixed_order
{
    //=============================================================================
    class PointQuadrature : public IntegrationRule
    {
    public:
        PointQuadrature()
            : IntegrationRule(0)
        {
        }

        IntegrationPoint point(int i) const override
        {
            (void)i;
            return {};
        }
    };
    //=============================================================================
    Point::Point()
        : NodalFiniteElement(mesh::CellType::point, 1,
                             std::make_unique<PointQuadrature>())
    {
    }
    //=============================================================================
    void Point::eval_shape(const std::array<real_t, 3> &pt,
                           la::DenseMatrix &N) const
    {
        (void)pt;
        N(0, 0) = 1.0;
    }
    //=============================================================================
    void Point::eval_shape_grad(const std::array<real_t, 3> &pt,
                                la::DenseMatrix &dNdxi) const
    {
        (void)pt;
        (void)dNdxi;
    }
}