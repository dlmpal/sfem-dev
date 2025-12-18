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

        FEData transform(int elem_idx, int pdim, const std::array<real_t, 3> &pt,
                         std::span<const std::array<real_t, 3>> elem_pts) const override;
    };
}