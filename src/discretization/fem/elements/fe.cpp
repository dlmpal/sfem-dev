#include "fe.hpp"

namespace sfem::fem
{
    //=============================================================================
    FiniteElement::FiniteElement(mesh::CellType cell_type, int order,
                                 std::unique_ptr<IntegrationRule> &&integration_rule)
        : cell_type_(cell_type),
          order_(order),
          integration_rule_(std::move(integration_rule))
    {
        if (order_ <= 0)
        {
            SFEM_ERROR(std::format("Invalid order {} (<=0)\n", order_));
        }
    }
    //=============================================================================
    mesh::CellType FiniteElement::cell_type() const
    {
        return cell_type_;
    }
    //=============================================================================
    int FiniteElement::order() const
    {
        return order_;
    }
    //=============================================================================
    IntegrationRule *FiniteElement::integration_rule() const
    {
        return integration_rule_.get();
    }
    //=============================================================================
    int FiniteElement::n_nodes() const
    {
        return mesh::cell_num_nodes(cell_type_);
    }
    //=============================================================================
    int FiniteElement::dim() const
    {
        return mesh::cell_dim(cell_type_);
    }
    //=============================================================================
    FEData::FEData(int n_nodes, int pdim, int gdim)
        : detJ(0.0),
          dXdxi(pdim, gdim),
          dxidX(gdim, pdim),
          N(n_nodes, 1),
          dNdxi(n_nodes, gdim),
          dNdX(n_nodes, pdim)
    {
    }
}
