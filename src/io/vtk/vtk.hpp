#pragma once

#include "../../mesh/mesh.hpp"
#include "../../discretization/fem/fe_function.hpp"
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
    /// @param cell_funcs Cell function data
    /// @param node_funcs Nodal function data (i.e. point data)
    /// @param type VTK file type (legacy for .vtk and xml for .vtu)
    void write(std::filesystem::path filename,
               const std::vector<int> &cell_types,
               const graph::Connectivity &cell_to_node,
               const std::vector<std::array<real_t, 3>> &points,
               const std::vector<std::shared_ptr<const Function>> &cell_funcs,
               const std::vector<std::shared_ptr<const Function>> &node_funcs,
               VTKFileType type = VTKFileType::xml);

    /// @brief Export a mesh to VTK
    /// @param filename Output filename
    /// @param mesh Mesh
    /// @param cell_funcs Cell function data
    /// @param node_funcs Nodal function data (i.e. point data)
    /// @param type VTK file type
    void write(std::filesystem::path filename,
               const mesh::Mesh &mesh,
               const std::vector<std::shared_ptr<const Function>> &cell_funcs = {},
               const std::vector<std::shared_ptr<const Function>> &node_funcs = {},
               VTKFileType type = VTKFileType::xml);

    /// @brief Export a finite element function to VTK
    /// @param filename Output filename
    /// @param funcs Finite element function(s)
    /// @param type VTK file type
    /// @note All functions should belong to the same FE space
    void write(const std::filesystem::path &filename,
               const std::vector<std::shared_ptr<const fem::FEFunction>> &funcs,
               VTKFileType type = VTKFileType::xml);
}