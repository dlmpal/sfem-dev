#pragma once

#include "../../mesh/mesh.hpp"
#include "../../discretization/fem/fe_space.hpp"
#include <filesystem>

namespace sfem::io::vtk
{
    /// @brief Available VTK file types.
    /// All available types refer to unstructured grid formats
    enum class VTKFileType
    {
        legacy,
        xml
    };

    /// @brief Create either a modern (.vtu) or legacy (.vtk) VTK file
    /// @param filename Output filename
    /// @param cells Cells
    /// @param cell_orders Order of each cell
    /// @param cell_to_node Cell-to-node connectivity
    /// @param points xyz coordinates of each node
    /// @param cell_names Name for each set of cell values
    /// @param cell_values Cell values
    /// @param cell_names Name for each set of node values
    /// @param node_values Node values
    /// @param type VTK file type (legacy for .vtk and xml for .vtu)
    /// @todo Document filename behaviour
    void write(std::filesystem::path filename,
               const std::vector<mesh::Cell> &cells,
               const std::vector<int> &cell_orders,
               const graph::Connectivity &cell_to_node,
               const std::vector<std::array<real_t, 3>> &points,
               const std::vector<std::vector<real_t>> &cell_values = {},
               const std::vector<std::string> &cell_names = {},
               const std::vector<std::vector<real_t>> &node_values = {},
               const std::vector<std::string> &node_names = {},
               VTKFileType type = VTKFileType::xml);

    /// @brief Export a mesh to VTK
    /// @param filename Output filename
    /// @param mesh Mesh
    /// @param type VTK file type
    void write(std::filesystem::path filename,
               const mesh::Mesh &mesh,
               VTKFileType type = VTKFileType::xml);

    /// @brief Export a finite element space (and its values) to VTK
    /// @param filename Output filename
    /// @param fe_space Finite element space
    /// @param values Function values
    /// @param type VTK file type
    void write(const std::filesystem::path &filename,
               const fem::FESpace &fe_space,
               const std::vector<real_t> &values,
               VTKFileType type = VTKFileType::xml);
}