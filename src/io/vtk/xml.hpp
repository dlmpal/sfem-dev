#pragma once

#include "../../mesh/mesh.hpp"
#include "../../discretization/fem/fe_space.hpp"
#include <filesystem>

/// @brief VTK XML file writers
namespace sfem::io::vtk::xml
{
    /// @brief Create a VTK unstructured grid file (.vtu)
    /// @param filename Output filename
    /// @param cells Cells
    /// @param cell_orders Order of each cell
    /// @param cell_to_node Cell-to-node connectivity
    /// @param points xyz coordinates of each node
    /// @param cell_names Name for each set of cell values
    /// @param cell_values Cell values
    /// @param cell_names Name for each set of node values
    /// @param node_values Node values
    void write_vtu(const std::filesystem::path &filename,
                   const std::vector<mesh::Cell> &cells,
                   const std::vector<int> &cell_orders,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<std::vector<real_t>> &cell_values = {},
                   const std::vector<std::string> &cell_names = {},
                   const std::vector<std::vector<real_t>> &node_values = {},
                   const std::vector<std::string> &node_names = {});

    /// @brief Create a parallel VTK unstructured grid file (.pvtu)
    void write_pvtu(const std::filesystem::path &filename,
                    const std::vector<std::filesystem::path> &sources,
                    const std::vector<std::string> &cell_names = {},
                    const std::vector<std::string> &node_names = {});
}