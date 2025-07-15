#pragma once

#include "../../mesh/mesh.hpp"
#include "../../discretization/function.hpp"
#include <filesystem>

/// @brief VTK legacy file writers
namespace sfem::io::vtk::legacy
{
    /// @brief Create a legacy VTK file (.vtk)
    /// @param filename Output filename
    /// @param cell_types VTK cell types
    /// @param cell_to_node Cell-to-node (or cell-to-point) connectivity
    /// @param points Point xyz coordinates
    /// @param cell_funcs Cell function data
    /// @param node_funcs Nodal function data (i.e. point data)
    /// @note Not MPI collective
    void write_vtk(const std::filesystem::path &filename,
                   const std::vector<int> &cell_types,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<std::shared_ptr<const Function>> &cell_funcs,
                   const std::vector<std::shared_ptr<const Function>> &node_funcs);
}