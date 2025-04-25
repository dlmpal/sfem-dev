#pragma once

#include "../../mesh/mesh.hpp"
#include "../../discretization/fem/fe_space.hpp"
#include <filesystem>

/// @brief VTK legacy file writers
namespace sfem::io::vtk::legacy
{
    /// @brief Create a legacy VTK file (.vtk)
    /// @param filename Output filename
    /// @param cell_types VTK cell types
    /// @param cell_to_node Cell-to-node connectivity
    /// @param points Point xyz coordinates
    /// @param cell_names Name for each set of cell values
    /// @param cell_values Cell values
    /// @param cell_names Name for each set of node values
    /// @param node_values Node values
    /// @note Not MPI collective
    void write_vtk(const std::filesystem::path &filename,
                   const std::vector<int> &cell_types,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<std::vector<real_t>> &cell_values = {},
                   const std::vector<std::string> &cell_names = {},
                   const std::vector<std::vector<real_t>> &node_values = {},
                   const std::vector<std::string> &node_names = {});
}