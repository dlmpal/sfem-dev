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
    /// @param cell_to_node Cell-to-node (or cell-to-point) connectivity
    /// @param points Point xyz coordinates
    /// @param cell_data Cell data
    /// @param node_data Node data (i.e. point data)
    /// @note Not MPI collective
    void write_vtk(const std::filesystem::path &filename,
                   const std::vector<int> &cell_types,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<std::pair<std::string, std::span<real_t>>> &cell_data,
                   const std::vector<std::pair<std::string, std::span<real_t>>> &node_data);
}