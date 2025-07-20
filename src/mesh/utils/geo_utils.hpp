#pragma once

#include "../cell.hpp"
#include "../../geo/vec3.hpp"

namespace sfem::mesh
{
    /// @brief Compute the facet tangent vectors
    /// @note If the facet is a node or an edge, returns the same vector twice
    std::array<geo::Vec3, 2> facet_tangents(CellType facet_type,
                                            std::span<const std::array<real_t, 3>> facet_points);

    /// @brief Compute the facet normal vector
    geo::Vec3 facet_normal(CellType facet_type,
                           std::span<const std::array<real_t, 3>> facet_points);

    /// @brief Compute the midpoint of a cell
    std::array<real_t, 3> cell_midpoint(std::span<const std::array<real_t, 3>> cell_points);
}