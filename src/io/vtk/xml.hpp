#pragma once

#include "../../mesh/mesh.hpp"
#include "../../discretization/function.hpp"
#include <filesystem>

/// @brief VTK XML file writers
namespace sfem::io::vtk::xml
{
    /// @brief Create a VTK unstructured grid file (.vtu)
    /// @param filename Output filename
    /// @param cell_types VTK cell types
    /// @param cell_to_node Cell-to-node (or cell-to-point) connectivity
    /// @param points Point xyz coordinates
    /// @param cell_funcs Cell function data
    /// @param node_funcs Nodal function data (i.e. point data)
    void write_vtu(const std::filesystem::path &filename,
                   const std::vector<int> &cell_types,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<std::shared_ptr<const Function>> &cell_funcs,
                   const std::vector<std::shared_ptr<const Function>> &node_funcs);

    /// @brief Create a parallel VTK unstructured grid file (.pvtu)
    void write_pvtu(const std::filesystem::path &filename,
                    const std::vector<std::filesystem::path> &sources,
                    const std::vector<std::shared_ptr<const Function>> &cell_funcs,
                    const std::vector<std::shared_ptr<const Function>> &node_funcs);
}