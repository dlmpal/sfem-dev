#pragma once

#include <sfem/mesh/cell.hpp>

namespace sfem::io::vtk
{
    /// @brief Get the VTK equivalent for a native cell type of given order
    int cell_type_to_vtk(mesh::CellType cell_type, int order);

    /// @brief Get the corresponding node ordering for VTK cells.
    void cell_node_ordering_to_vtk(int vtk_type, std::span<int> nodes);
}