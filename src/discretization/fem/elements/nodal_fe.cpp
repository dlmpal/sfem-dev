#include "nodal_fe.hpp"

namespace sfem::fem
{
    //=============================================================================
    NodalFiniteElement::NodalFiniteElement(mesh::CellType cell_type, int order,
                                           std::unique_ptr<IntegrationRule> &&integration_rule)
        : FiniteElement(cell_type, order, std::move(integration_rule))
    {
    }
    //=============================================================================
    int NodalFiniteElement::n_nodes() const
    {
        return dof::cell_num_dof(cell_type_, order_);
    }
    //=============================================================================
    FEData NodalFiniteElement::transform(int elem_idx, int pdim, const std::array<real_t, 3> &pt,
                                         std::span<const std::array<real_t, 3>> elem_pts) const
    {
        FEData data(elem_idx, n_nodes(), pdim, dim());
        data.pt = pt;

        // Evaluate the shape function and its gradient w.r.t natural coordinates
        eval_shape(pt, data.N);
        eval_shape_grad(pt, data.dNdxi);

        // Return early for point elements
        if (dim() == 0)
        {
            data.detJ = 1.0;
            return data;
        }

        // Evaluate the natural to physical Jacobian
        for (int i = 0; i < pdim; i++)
        {
            for (int j = 0; j < dim(); j++)
            {
                for (int k = 0; k < n_nodes(); k++)
                {
                    data.dXdxi(i, j) += data.dNdxi(k, j) * elem_pts[k][i];
                }
            }
        }

        // Evaluate the physical to natural, or inverse, Jacobian
        // and the determinant
        std::tie(data.dxidX, data.detJ) = data.dXdxi.invert();

        // Check for non-positive jacobian
        if (data.detJ <= 0)
        {
            /// @todo
            SFEM_ERROR(std::format("Negative Jacobian"));
        }

        // Evaluate the shape function gradient w.r.t physical coordinates
        data.dNdX = data.dNdxi * data.dxidX;

        return data;
    }
}