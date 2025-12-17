#pragma once

#include "../io_field.hpp"
#include "../../mesh/mesh.hpp"
#include <filesystem>

/// @brief VTK legacy file writers
namespace sfem::io::vtk::legacy
{
    /// @brief Create a legacy VTK file (.vtk)
    /// @param filename Output filename
    /// @param cell_types VTK cell types
    /// @param cell_to_node Cell-to-node (or cell-to-point) connectivity
    /// @param points Point xyz coordinates
    /// @param cell_fields Cell field data
    /// @param node_fields Nodal field data (i.e. point data)
    /// @note Not MPI collective
    void write_vtk(const std::filesystem::path &filename,
                   const std::vector<int> &cell_types,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<IOField> &cell_fields,
                   const std::vector<IOField> &node_fields);
}