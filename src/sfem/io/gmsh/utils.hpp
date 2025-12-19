#pragma once

#include <sfem/mesh/cell.hpp>

namespace sfem::io::gmsh
{
    /// @brief Get the cell type corresponding to a Gmsh cell type
    /// @param gmsh_type Gmsh cell type
    /// @param gmsh_idx Gmsh index of the cell (used in error messages)
    /// @return Corresponding Cell type
    mesh::CellType gmsh_type_to_native(int gmsh_type, int gmsh_idx);
}