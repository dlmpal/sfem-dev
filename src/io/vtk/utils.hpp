#pragma once

#include "../../mesh/cell.hpp"

namespace sfem::io::vtk
{
    /// @brief Get the VTK equivalent for a native cell type of given order
    int cell_type_to_vtk(mesh::CellType cell_type, int order);
}