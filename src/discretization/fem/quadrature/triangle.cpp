#include "triangle.hpp"
#include "../../../base/error.hpp"

namespace sfem::fem::quadrature
{
    //=============================================================================
    static real_t triangle_qweights(int n_points, int i)
    {
        static const real_t triangle_qwt_1[] = {1.0};

        static const real_t triangle_qwt_3[] = {1.0 / 6.0,
                                                1.0 / 6.0,
                                                1.0 / 6.0};

        static const real_t triangle_qwt_4[] = {-27.0 / 96.0,
                                                25.0 / 96.0,
                                                25.0 / 96.0,
                                                25.0 / 96.0};

        static const real_t triangle_qwt_6[] = {0.054975871827661,
                                                0.054975871827661,
                                                0.054975871827661,
                                                0.111690794839006,
                                                0.111690794839006,
                                                0.111690794839006};
        switch (n_points)
        {
        case 1:
            return triangle_qwt_1[0];
        case 3:
            return triangle_qwt_3[i];
        case 4:
            return triangle_qwt_4[i];
        case 6:
            return triangle_qwt_6[i];
        default:
            SFEM_ERROR(std::format("Triangle integration not defined for {} points\n", n_points));
            return 0;
        }
    }
    //=============================================================================
    static std::array<real_t, 3>
    triangle_qpoints(int n_points, int i)
    {
        static const real_t triangle_qpt_1[] = {1.0 / 3.0, 1.0 / 3.0};

        static const real_t triangle_qpt_3[] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0,
                                                1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};

        static const real_t triangle_qpt_4[] = {1.0 / 3.0, 1.0 / 3.0,
                                                1.0 / 5.0, 1.0 / 5.0,
                                                3.0 / 5.0, 1.0 / 5.0,
                                                1.0 / 5.0, 3.0 / 5.0};

        static const real_t triangle_qpt_6[] = {
            0.091576213509771, 0.091576213509771,
            0.816847572980459, 0.091576213509771,
            0.091576213509771, 0.816847572980459,
            0.108103018168070, 0.108103018168070,
            0.445948490915965, 0.108103018168070,
            0.108103018168070, 0.445948490915965};

        switch (n_points)
        {
        case 1:
            return {triangle_qpt_1[0], triangle_qpt_1[1]};
        case 3:
            return {triangle_qpt_3[i * 2 + 0], triangle_qpt_3[i * 2 + 1], 0};
        case 4:
            return {triangle_qpt_4[i * 2 + 0], triangle_qpt_4[i * 2 + 1], 0};
        case 6:
            return {triangle_qpt_6[i * 2 + 0], triangle_qpt_6[i * 2 + 1], 0};
        default:
            SFEM_ERROR(std::format("Triangle integration not defined for {} points\n", n_points));
            return {};
        }
    }
    //=============================================================================
    Triangle::Triangle(int n_points)
        : IntegrationRule(n_points)
    {
    }
    //=============================================================================
    real_t Triangle::weight(int qpt_idx) const
    {
        return triangle_qweights(n_points_, qpt_idx);
    }
    //=============================================================================
    std::array<real_t, 3> Triangle::point(int qpt_idx) const
    {
        return triangle_qpoints(n_points_, qpt_idx);
    }
}