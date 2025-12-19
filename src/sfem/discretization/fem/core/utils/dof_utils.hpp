#pragma once

#include <sfem/mesh/cell.hpp>
#include <vector>

/// @brief Miscellaneous utilities related to DoF handling
namespace sfem::fem::dof
{
    /// @brief Get the number of DoF per cell type and polynomial order
    int cell_num_dof(mesh::CellType cell_type, int order);

    /// @brief Get the number of internal DoF per cell type and polynomial order
    int cell_num_internal_dof(mesh::CellType cell_type, int order);

    /// @brief Compute the physical coordinates of the DoFs for a line segment
    /// @param order The polynomial order
    /// @param points The line segment start and end points
    void compute_line_dof_points(int order, std::vector<std::array<real_t, 3>> &points);

    /// @brief Compute the physical coordinates of the DoFs for a triangle
    /// @param order The polynomial order
    /// @param points The triangle corner node coordinates
    void compute_triangle_dof_points(int order, std::vector<std::array<real_t, 3>> &points);

    /// @brief Compute the physical coordinates of the DoFs for a quadrilateral
    /// @param order The polynomial order
    /// @param points The quadrilateral corner node coordinates
    void compute_quad_dof_points(int order, std::vector<std::array<real_t, 3>> &points);

    /// @brief Compute the physical coordinates of the DoFs for a tetrahedron
    /// @param order The polynomial order
    /// @param points The tetrahedron corner node coordinates
    void compute_tet_dof_points(int order, std::vector<std::array<real_t, 3>> &points);

    /// @brief Compute the physical coordinates of the DoFs for a hexahedron
    /// @param order The polynomial order
    /// @param points The hexahedron corner node coordinates
    void compute_hex_dof_points(int order, std::vector<std::array<real_t, 3>> &points);

    /// @brief Compute the physical coordinates of the DoFs for a given
    /// cell and polynomial order. For example:
    //
    // Quadrilateral with order = 3
    // 3---9---8---2
    // |           |
    // 10  14  15  7
    // |           |
    // 11  12  13  6
    // |           |
    // 0---4---5---1
    //
    // Triangle with order = 3
    // 2
    // |  \
    // 7    6
    // |      \
    // 8   9   5
    // |         \
    // 0---3---4---1
    //
    /// @param cell_type The cell type
    /// @param order The polynomial order
    /// @param points the cell's corner node coordinates.
    /// For order > 1, this vector will be extended by the
    /// coordinates of the edge, face and internal DoFs
    ///
    /// @todo Add support for higher order elements in 3D
    void compute_cell_dof_points(mesh::CellType cell_type, int order,
                                 std::vector<std::array<real_t, 3>> &points);
}