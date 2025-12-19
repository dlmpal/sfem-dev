#pragma once

#include <sfem/geo/vec3.hpp>

namespace sfem::geo
{
    /// @brief Compute the distance between two points
    real_t compute_distance(const std::array<real_t, 3> &point1,
                            const std::array<real_t, 3> &point2);

    /// @brief Compute the midpoint of a line segment
    std::array<real_t, 3> compute_line_midpoint(const std::array<real_t, 3> &X1,
                                                const std::array<real_t, 3> &X2);

    /// @brief Compute the xyz coordinates of the n-th point, on an
    /// line segment divided into k equal sized parts
    std::array<real_t, 3> compute_line_nth_point(int n, int k,
                                                 const std::array<real_t, 3> &X1,
                                                 const std::array<real_t, 3> &X2);

    /// @brief Compute the intersection of two line segments
    /// @param X1 Start point for the first segment
    /// @param X2 End point for the first segment
    /// @param X3 Start point for the second segment
    /// @param X4 End point for the second segment
    /// @return The point of intersection
    std::array<real_t, 3> compute_line_intersection(const std::array<real_t, 3> &X1,
                                                    const std::array<real_t, 3> &X2,
                                                    const std::array<real_t, 3> &X3,
                                                    const std::array<real_t, 3> &X4);

    /// @brief Compute the angle between two vectors
    /// @param v1 First vector
    /// @param v2 Second vector
    /// @return Angle (in radians)
    real_t vector_angle(const Vec3 &v1, const Vec3 &v2);
}