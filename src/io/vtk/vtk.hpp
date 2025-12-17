#pragma once

#include "../io_field.hpp"
#include "../../mesh/mesh.hpp"
#include "../../discretization/fem/fe_field.hpp"
#include "../../discretization/fvm/fv_field.hpp"
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
    /// @param cell_fields Cell field data
    /// @param node_fields Nodal field data (i.e. point data)
    /// @param type VTK file type (legacy for .vtk and xml for .vtu)
    void write(std::filesystem::path filename,
               const std::vector<int> &cell_types,
               const graph::Connectivity &cell_to_node,
               const std::vector<std::array<real_t, 3>> &points,
               const std::vector<IOField> &cell_fields,
               const std::vector<IOField> &node_fields,
               VTKFileType type = VTKFileType::xml);

    /// @brief Export a mesh to VTK
    /// @param filename Output filename
    /// @param mesh Mesh
    /// @param cell_fields Cell field data
    /// @param node_fields Nodal field data (i.e. point data)
    /// @param type VTK file type
    void write(std::filesystem::path filename,
               const mesh::Mesh &mesh,
               const std::vector<IOField> &cell_fields = {},
               const std::vector<IOField> &node_fields = {},
               VTKFileType type = VTKFileType::xml);

    /// @brief Export a set of finite volume fields to VTK
    /// @param filename Output filename
    /// @param fields Finite volume fields
    /// @param type VTK file type
    void write(const std::filesystem::path &filename,
               const std::vector<fvm::FVField> &fields,
               VTKFileType type = VTKFileType::xml);

    /// @brief Export a set of finite element field to VTK
    /// @param filename Output filename
    /// @param fields Finite element fields
    /// @param type VTK file type
    /// @note All fields should belong to the same FE space
    void write(const std::filesystem::path &filename,
               const std::vector<fem::FEField> &nodal_fields,
               const std::vector<fem::FEField> &cell_fields = {},
               VTKFileType type = VTKFileType::xml);
}