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
}