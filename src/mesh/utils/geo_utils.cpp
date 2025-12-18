#include "geo_utils.hpp"

namespace sfem::mesh
{
    //=============================================================================
    std::array<geo::Vec3, 2>
    facet_tangents(CellType facet_type,
                   std::span<const std::array<real_t, 3>> facet_points)
    {
        if (facet_type == CellType::point)
        {
            return {};
        }
        else if (facet_type == CellType::line)
        {
            return {geo::Vec3(facet_points[0], facet_points[1]),
                    geo::Vec3(facet_points[0], facet_points[1])};
        }
        else
        {
            return {geo::Vec3(facet_points[0], facet_points[1]),
                    geo::Vec3(facet_points[0], facet_points[2])};
        }
    }
    //=============================================================================
    geo::Vec3 facet_normal(CellType facet_type,
                           std::span<const std::array<real_t, 3>> facet_points)
    {
        if (facet_type == CellType::point)
        {
            return geo::Vec3(1, 0, 0);
        }
        else if (facet_type == CellType::line)
        {
            auto [t1, _] = facet_tangents(facet_type, facet_points);
            return t1.normal2D();
        }
        else
        {
            auto [t1, t2] = facet_tangents(facet_type, facet_points);
            return geo::cross(t1, t2);
        }
    }
    //=============================================================================
    std::array<real_t, 3> cell_midpoint(std::span<const std::array<real_t, 3>> cell_points)
    {
        real_t w = 1.0 / static_cast<real_t>(cell_points.size());
        std::array<real_t, 3> midpoint{};
        for (const auto &point : cell_points)
        {
            for (int i = 0; i < 3; i++)
            {
                midpoint[i] += point[i] * w;
            }
        }
        return midpoint;
    }
    //=============================================================================
    std::array<real_t, 3> map_facet_to_cell_ref(CellType cell_type, int facet_idx,
                                                const std::array<real_t, 3> &xi_facet)
    {
        const real_t xi = xi_facet[0];
        const real_t eta = xi_facet[1];

        switch (cell_type)
        {
        case CellType::tetrahedron:
        {
            switch (facet_idx)
            {
            case 0:
                return {xi, eta, 1.0 - xi - eta};
            case 1:
                return {0.0, xi, eta};
            case 2:
                return {xi, 0.0, eta};
            case 3:
                return {xi, eta, 0.0};
            default:
                SFEM_ERROR(std::format("Invalid facet index ({}) for tetrahedron\n", facet_idx));
                return {};
            }
        }
        case CellType::hexahedron:
        {
            switch (facet_idx)
            {
            case 0:
                return {-1.0, xi, eta};
            case 1:
                return {1.0, xi, eta};
            case 2:
                return {xi, -1.0, eta};
            case 3:
                return {xi, 1.0, eta};
            case 4:
                return {xi, eta, -1.0};
            case 5:
                return {xi, eta, 1.0};
            default:
                SFEM_ERROR(std::format("Invalid facet index ({}) for hexahedron\n", facet_idx));
                return {};
            }
        }
        case CellType::triangle:
        {
            switch (facet_idx)
            {
            case 0:
                return {xi, 0.0, 0.0};
            case 1:
                return {1.0 - xi, xi, 0.0};
            case 2:
                return {0.0, 1.0 - xi, 0.0};
            default:
                SFEM_ERROR(std::format("Invalid facet index ({}) for triangle\n", facet_idx));
                return {};
            }
        }
        case CellType::quadrilateral:
        {
            switch (facet_idx)
            {
            case 0:
                return {xi, -1.0, 0.0};
            case 1:
                return {1.0, xi, 0.0};
            case 2:
                return {xi, 1.0, 0.0};
            case 3:
                return {-1.0, xi, 0.0};
            default:
                SFEM_ERROR(std::format("Invalid facet index ({}) for quadrilateral\n", facet_idx));
                return {};
            }
        }
        default:
            SFEM_BAD_CELL_ERROR(cell_type);
            return {};
        }
    }
}