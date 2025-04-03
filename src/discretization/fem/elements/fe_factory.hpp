#pragma once

#include "fixed_order.hpp"

namespace sfem::fem
{
    /// @brief Create a nodal finite element for a given cell type and order
    std::shared_ptr<NodalFiniteElement>
    create_nodal_element(mesh::CellType cell_type, int order);
}