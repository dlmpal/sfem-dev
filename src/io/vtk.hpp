#pragma once

#include "../mesh/mesh.hpp"
#include "../discretization/fem/fe_space.hpp"
#include <filesystem>

/// @brief VTK file I/O
namespace sfem::io::vtk
{
    /// @brief Legacy VTK file writer
    /// @param filename Output filename
    /// @param cells The cells
    /// @param cell_orders The order of each cell
    /// @param cell_to_node The cell-to-node connectivity
    /// @param points The xyz coordinates of each node
    /// @param cell_names The name for each set of cell values
    /// @param cell_values Cell values
    /// @param cell_names The name for each set of node values
    /// @param node_values Node values
    /// @note Not MPI collective
    void write(const std::filesystem::path &filename,
               const std::vector<mesh::Cell> &cells,
               const std::vector<int> &cell_orders,
               const graph::Connectivity &cell_to_node,
               const std::vector<std::array<real_t, 3>> &points,
               const std::vector<std::vector<real_t>> &cell_values = {},
               const std::vector<std::string> &cell_names = {},
               const std::vector<std::vector<real_t>> &node_values = {},
               const std::vector<std::string> &node_names = {});

    /// @brief Export mesh to VTK
    void write(const std::filesystem::path &filename,
               const mesh::Mesh &mesh);

    void write(const std::filesystem::path &filename,
               const fem::FESpace &fe_space,
               const std::vector<real_t> &values);
}