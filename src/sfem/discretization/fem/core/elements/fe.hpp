#pragma once

#include <sfem/discretization/fem/core/quadrature/quadrature.hpp>
#include <sfem/discretization/fem/core/utils/dof_utils.hpp>
#include <sfem/la/native/dense_matrix.hpp>
#include <memory>
#include <functional>

namespace sfem::fem
{
    // Forward declaration
    class FEData;

    /// @brief Finite element ABC
    class FiniteElement
    {
    public:
        FiniteElement(mesh::CellType cell_type, int order,
                      std::unique_ptr<IntegrationRule> &&integration_rule);

        virtual ~FiniteElement() = default;

        /// @brief Get the element's reference cell type
        mesh::CellType cell_type() const;

        /// @brief Get the element's order (i.e. polynomial degree)
        int order() const;

        /// @brief Get the integration rule
        IntegrationRule *integration_rule() const;

        /// @brief Get the number of nodes (i.e. DoF) for the element
        virtual int n_nodes() const;

        /// @brief Get the topological dimension of the element's reference cell type
        int dim() const;

        /// @brief Evaluate the element coordinate transform for a given point
        /// @note The physical dimension must be greater than or equal to the element's reference
        /// dimension
        virtual FEData transform(int elem_idx, int pdim, const std::array<real_t, 3> &pt,
                                 std::span<const std::array<real_t, 3>> elem_pts) const = 0;

        /// @brief Evaluate the element's shape functions at a given point
        virtual void eval_shape(const std::array<real_t, 3> &pt,
                                la::DenseMatrix &N) const = 0;

        /// @brief Evaluate the gradients of the element's shape functions at a given point
        virtual void eval_shape_grad(const std::array<real_t, 3> &pt,
                                     la::DenseMatrix &dNdx) const = 0;

    protected:
        /// @brief The element's reference cell type
        mesh::CellType cell_type_;

        /// @brief The element's order (i.e. polynomial degree)
        int order_;

        /// @brief The integration rule (or quadrature)
        std::unique_ptr<IntegrationRule> integration_rule_;
    };

    /// @brief Finite element coordinate transform data
    struct FEData
    {
        FEData(int elem_idx, int n_nodes, int pdim, int gdim);

        /// @brief Current element index
        int elem_idx;

        /// @brief Number of nodes
        int n_nodes;

        /// @brief Physical dimension
        int pdim;

        /// @brief Geometric/topological dimension
        int gdim;

        /// @brief Transform evaluation point
        std::array<real_t, 3> pt;

        /// @brief Natural-to-physical Jacobian determinant
        real_t detJ;

        /// @brief Natural-to-physical Jacobian (direct)
        la::DenseMatrix dXdxi;

        /// @brief Physical-to-natural jacobian (inverse)
        la::DenseMatrix dxidX;

        /// @brief Shape function
        la::DenseMatrix N;

        /// @brief Shape function gradient (natural)
        la::DenseMatrix dNdxi;

        /// @brief Shape function gradient (physical)
        la::DenseMatrix dNdX;
    };

    using ElementOperator = std::function<void(const FEData &data, la::DenseMatrix &f)>;
}