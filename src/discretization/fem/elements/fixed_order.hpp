#pragma once

#include "fe.hpp"
#include "../quadrature/gauss.hpp"
#include "../quadrature/triangle.hpp"
#include "../quadrature/tetrahedron.hpp"

/// @brief Fixed-order finite elements
namespace sfem::fem::fixed_order
{
    class Point : public NodalFiniteElement
    {
    public:
        class PointQuadrature : public IntegrationRule
        {
        public:
            PointQuadrature()
                : IntegrationRule(0)
            {
            }

            int n_points() const override
            {
                return 1;
            }

            IntegrationPoint point(int i) const override
            {
                (void)i;
                return {};
            }
        };

        Point()
            : NodalFiniteElement(mesh::CellType::point, 1,
                                 std::make_unique<PointQuadrature>())
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override
        {
            (void)pt;
            N(0, 0) = 1.0;
        }

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override
        {
            (void)pt;
            (void)dNdxi;
        }
    };

    class Line2 : public NodalFiniteElement
    {
    public:
        Line2()
            : NodalFiniteElement(mesh::CellType::line, 1,
                                 std::make_unique<quadrature::Gauss1D>(1))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Line3 : public NodalFiniteElement
    {
    public:
        Line3()
            : NodalFiniteElement(mesh::CellType::line, 2,
                                 std::make_unique<quadrature::Gauss1D>(2))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Line4 : public NodalFiniteElement
    {
    public:
        Line4()
            : NodalFiniteElement(mesh::CellType::line, 3,
                                 std::make_unique<quadrature::Gauss1D>(3))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tri3 : public NodalFiniteElement
    {
    public:
        Tri3()
            : NodalFiniteElement(mesh::CellType::triangle, 1,
                                 std::make_unique<quadrature::Triangle>(1))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tri6 : public NodalFiniteElement
    {
    public:
        Tri6()
            : NodalFiniteElement(mesh::CellType::triangle, 2,
                                 std::make_unique<quadrature::Triangle>(2))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tri10 : public NodalFiniteElement
    {
    public:
        Tri10()
            : NodalFiniteElement(mesh::CellType::triangle, 3,
                                 std::make_unique<quadrature::Triangle>(3))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Quad4 : public NodalFiniteElement
    {
    public:
        Quad4()
            : NodalFiniteElement(mesh::CellType::quad, 1,
                                 std::make_unique<quadrature::Gauss2D>(1))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Quad9 : public NodalFiniteElement
    {
    public:
        Quad9()
            : NodalFiniteElement(mesh::CellType::quad, 2,
                                 std::make_unique<quadrature::Gauss2D>(2))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Quad16 : public NodalFiniteElement
    {
    public:
        Quad16()
            : NodalFiniteElement(mesh::CellType::quad, 3,
                                 std::make_unique<quadrature::Gauss2D>(3))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tet4 : public NodalFiniteElement
    {
    public:
        Tet4()
            : NodalFiniteElement(mesh::CellType::tet, 1,
                                 std::make_unique<quadrature::Tetrahedron>(1))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tet10 : public NodalFiniteElement
    {
    public:
        Tet10()
            : NodalFiniteElement(mesh::CellType::tet, 2,
                                 std::make_unique<quadrature::Tetrahedron>(2))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Hex8 : public NodalFiniteElement
    {
    public:
        Hex8()
            : NodalFiniteElement(mesh::CellType::hex, 1,
                                 std::make_unique<quadrature::Gauss3D>(1))
        {
        }

        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;

        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };
}