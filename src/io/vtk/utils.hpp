#pragma once

#include "../../mesh/cell.hpp"

namespace sfem::io::vtk
{
    int cell_type_to_vtk(mesh::CellType cell_type, int order);
}