#include "utils.hpp"
#include "../la/native/dense_matrix_utils.hpp"
#include <cmath>

namespace sfem::geo
{
    //=============================================================================
    real_t compute_distance(const std::array<real_t, 3> &point1,
                            const std::array<real_t, 3> &point2)
    {
        real_t dx = point2[0] - point1[0];
        real_t dy = point2[1] - point1[1];
        real_t dz = point2[2] - point1[2];
        return std::sqrt(dx * dx +
                         dy * dy +
                         dz * dz);
    }
    //=============================================================================
    std::array<real_t, 3> compute_line_midpoint(const std::array<real_t, 3> &X1,
                                                const std::array<real_t, 3> &X2)
    {
        return {(X1[0] + X2[0]) * 0.5,
                (X1[1] + X2[1]) * 0.5,
                (X1[2] + X2[2]) * 0.5};
    }
    //=============================================================================
    std::array<real_t, 3> compute_line_nth_point(int n, int k,
                                                 const std::array<real_t, 3> &X1,
                                                 const std::array<real_t, 3> &X2)
    {
        real_t dx = (X2[0] - X1[0]) / k;
        real_t dy = (X2[1] - X1[1]) / k;
        real_t dz = (X2[2] - X1[2]) / k;

        return {X1[0] + dx * n,
                X1[1] + dy * n,
                X1[2] + dz * n};
    }
    //=============================================================================
    std::array<real_t, 3> compute_line_intersection(const std::array<real_t, 3> &X1,
                                                    const std::array<real_t, 3> &X2,
                                                    const std::array<real_t, 3> &X3,
                                                    const std::array<real_t, 3> &X4)
    {
        real_t dX1 = X2[0] - X1[0];
        real_t dY1 = X2[1] - X1[1];
        real_t dZ1 = X2[2] - X1[2];

        real_t dX2 = X4[0] - X3[0];
        real_t dY2 = X4[1] - X3[1];
        real_t dZ2 = X4[2] - X3[2];

        std::array<real_t, 3 * 2> A = {dX1, -dX2,
                                       dY1, -dY2,
                                       dZ1, -dZ2};

        std::array<real_t, 3> b = {X3[0] - X1[0],
                                   X3[1] - X1[1],
                                   X3[2] - X1[2]};

        std::array<real_t, 2 * 3> Ainv;
        la::utils::pinv(3, 2, A, Ainv);

        std::array<real_t, 2> x;
        la::utils::matmult(2, 1, 3, Ainv, b, x);

        return {X1[0] + dX1 * x[0],
                X1[1] + dY1 * x[0],
                X1[2] + dZ1 * x[0]};
    }
    //=============================================================================
    real_t vector_angle(const Vec3 &v1, const Vec3 &v2)
    {
        real_t mag1 = v1.mag();
        real_t mag2 = v2.mag();
        real_t prod = inner(v1, v2);
        return std::acos(prod / mag1 / mag2);
    }
}