#pragma once

#include "../quadrature/quadrature.hpp"
#include "../utils/dof_utils.hpp"
#include "../../../la/native/dense_matrix.hpp"
#include <memory>

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

        /// @brief Compute the coordinate transform for the element, for a given
        /// physical dimension
        /// @note The physical dimension must be greater than or equal to the element's reference
        /// dimension
        virtual FEData transform(int pdim, const IntegrationPoint &qpoint,
                                 std::span<const std::array<real_t, 3>> points) const = 0;

        /// @brief Evaluate the element's shape functions at a given point
        virtual void eval_shape(const std::array<real_t, 3> &point,
                                la::DenseMatrix &N) const = 0;

        /// @brief Evaluate the gradients of the element's shape functions at a given point
        virtual void eval_shape_grad(const std::array<real_t, 3> &point,
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
        FEData(int n_nodes, int pdim, int gdim);

        /// @brief Natural-to-physical Jacobian determinant
        /// @note Multiplied by the quadrature weight
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
}