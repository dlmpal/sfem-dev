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
        FiniteElement(mesh::CellType cell_type,
                      int order,
                      std::unique_ptr<IntegrationRule> &&integration_rule)
            : cell_type_(cell_type),
              order_(order),
              integration_rule_(std::move(integration_rule))
        {
            if (order_ <= 0)
            {
                /// @todo error
            }
        }

        virtual ~FiniteElement() = default;

        mesh::CellType cell_type() const
        {
            return cell_type_;
        }

        int order() const
        {
            return order_;
        }

        IntegrationRule *integration_rule() const
        {
            return integration_rule_.get();
        }

        virtual int n_nodes() const
        {
            return mesh::cell_num_nodes(cell_type_);
        }

        int dim() const
        {
            return mesh::cell_dim(cell_type_);
        }

        virtual FEData transform(int pdim,
                                 const IntegrationPoint &qpoint,
                                 std::span<const std::array<real_t, 3>> points) const = 0;

        virtual void eval_shape(const std::array<real_t, 3> &pt,
                                la::DenseMatrix &N) const = 0;

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

    /// @brief The (classic) finite element where the DoF are located
    /// at the nodes of the cell
    class NodalFiniteElement : public FiniteElement
    {
    public:
        NodalFiniteElement(mesh::CellType cell_type,
                           int order,
                           std::unique_ptr<IntegrationRule> &&integration_rule)
            : FiniteElement(cell_type, order, std::move(integration_rule))
        {
        }

        int n_nodes() const override
        {
            return dof::cell_num_dof(cell_type_, order_);
        }

        FEData transform(int pdim,
                         const IntegrationPoint &qpoint,
                         std::span<const std::array<real_t, 3>> points) const override;
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