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
    /// @param cell_types VTK cell types
    /// @param cell_to_node Cell-to-node connectivity
    /// @param points Point xyz coordinates
    /// @param cell_data Cell data
    /// @param node_data Node data (i.e. point data)
    /// @param type VTK file type (legacy for .vtk and xml for .vtu)
    /// @todo Document filename behaviour
    void write(std::filesystem::path filename,
               const std::vector<int> &cell_types,
               const graph::Connectivity &cell_to_node,
               const std::vector<std::array<real_t, 3>> &points,
               const std::vector<std::pair<std::string, std::span<real_t>>> &cell_data,
               const std::vector<std::pair<std::string, std::span<real_t>>> &node_data,
               VTKFileType type = VTKFileType::xml);

    /// @brief Export a mesh to VTK
    /// @param filename Output filename
    /// @param mesh Mesh
    /// @param cell_data Cell data
    /// @param node_data Node data (i.e. point data)
    /// @param type VTK file type
    void write(std::filesystem::path filename,
               const mesh::Mesh &mesh,
               const std::vector<std::pair<std::string, std::span<real_t>>> &cell_data = {},
               const std::vector<std::pair<std::string, std::span<real_t>>> &node_data = {},
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