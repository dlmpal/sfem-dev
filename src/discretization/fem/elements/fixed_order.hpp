#pragma once

#include "nodal_fe.hpp"
#include "../quadrature/gauss.hpp"
#include "../quadrature/triangle.hpp"
#include "../quadrature/tetrahedron.hpp"

/// @brief Fixed-order (nodal) finite elements
namespace sfem::fem::fixed_order
{
    class Point : public NodalFiniteElement
    {
    public:
        Point();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Line2 : public NodalFiniteElement
    {
    public:
        Line2();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Line3 : public NodalFiniteElement
    {
    public:
        Line3();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Line4 : public NodalFiniteElement
    {
    public:
        Line4();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tri3 : public NodalFiniteElement
    {
    public:
        Tri3();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tri6 : public NodalFiniteElement
    {
    public:
        Tri6();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tri10 : public NodalFiniteElement
    {
    public:
        Tri10();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Quad4 : public NodalFiniteElement
    {
    public:
        Quad4();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Quad9 : public NodalFiniteElement
    {
    public:
        Quad9();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Quad16 : public NodalFiniteElement
    {
    public:
        Quad16();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tet4 : public NodalFiniteElement
    {
    public:
        Tet4();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Tet10 : public NodalFiniteElement
    {
    public:
        Tet10();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };

    class Hex8 : public NodalFiniteElement
    {
    public:
        Hex8();
        void eval_shape(const std::array<real_t, 3> &pt,
                        la::DenseMatrix &N) const override;
        void eval_shape_grad(const std::array<real_t, 3> &pt,
                             la::DenseMatrix &dNdxi) const override;
    };
}