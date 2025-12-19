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

        real_t weight(int) const override
        {
            return 1.0;
        }

        std::array<real_t, 3> point(int) const override
        {
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
    void Point::eval_shape(const std::array<real_t, 3> &,
                           la::DenseMatrix &N) const
    {
        N(0, 0) = 1.0;
    }
    //=============================================================================
    void Point::eval_shape_grad(const std::array<real_t, 3> &,
                                la::DenseMatrix &) const
    {
    }
}