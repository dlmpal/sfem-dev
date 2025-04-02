#pragma once

#include "fixed_order.hpp"

namespace sfem::fem
{
    std::shared_ptr<NodalFiniteElement>
    create_nodal_element(mesh::CellType cell_type, int order);
}