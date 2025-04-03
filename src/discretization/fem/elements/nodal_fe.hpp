#pragma once

#include "fe.hpp"

namespace sfem::fem
{
    /// @brief The (classic) finite element where the DoF are located
    /// at the nodes of the cell
    class NodalFiniteElement : public FiniteElement
    {
    public:
        NodalFiniteElement(mesh::CellType cell_type, int order,
                           std::unique_ptr<IntegrationRule> &&integration_rule);

        int n_nodes() const override;

        FEData transform(int pdim, const IntegrationPoint &qpoint,
                         std::span<const std::array<real_t, 3>> points) const override;
    };
}