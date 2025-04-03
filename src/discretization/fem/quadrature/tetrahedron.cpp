#include "tetrahedron.hpp"
#include "../../../base/error.hpp"

namespace sfem::fem::quadrature
{
    //=============================================================================
    static real_t tet_qweights(int n_points, int i)
    {
        static const real_t tet_qwt_4[] = {1.0 / 24.0,
                                           1.0 / 24.0,
                                           1.0 / 24.0,
                                           1.0 / 24.0};

        static const real_t tet_qwt_5[] = {-2.0 / 15.0,
                                           3.0 / 40.0,
                                           3.0 / 40.0,
                                           3.0 / 40.0,
                                           3.0 / 40.0};

        switch (n_points)
        {
        case 4:
            return tet_qwt_4[i];
        case 5:
            return tet_qwt_5[i];
        default:
            SFEM_ERROR(std::format("Tetrahedron integration not defined for {} points\n", n_points));
            return 0;
        }
    }
    //=============================================================================
    static std::array<real_t, 3>
    tet_qpoints(int n_points, int i)
    {
        static const real_t tet_qpt_4[] = {
            0.138196601125011, 0.138196601125011, 0.138196601125011,
            0.585410196624969, 0.138196601125011, 0.138196601125011,
            0.138196601125011, 0.585410196624969, 0.138196601125011,
            0.138196601125011, 0.138196601125011, 0.585410196624969};

        static const real_t tet_qpt_5[] = {0.25, 0.25, 0.25,
                                           1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0,
                                           0.5, 1.0 / 6.0, 1.0 / 6.0,
                                           1.0 / 6.0, 0.5, 1.0 / 6.0,
                                           1.0 / 6.0, 1.0 / 6.0, 0.5};
        switch (n_points)
        {
        case 4:
            return {tet_qpt_4[i * 3 + 0],
                    tet_qpt_4[i * 3 + 1],
                    tet_qpt_4[i * 3 + 2]};
        case 5:
            return {tet_qpt_5[i * 3 + 0],
                    tet_qpt_5[i * 3 + 1],
                    tet_qpt_5[i * 3 + 2]};
        default:
            SFEM_ERROR(std::format("Tetrahedron integration not defined for {} points\n", n_points));
            return {};
        }
    }
    //=============================================================================
    Tetrahedron::Tetrahedron(int n_points)
        : IntegrationRule(n_points)
    {
    }
    //=============================================================================
    IntegrationPoint Tetrahedron::point(int i) const
    {
        return IntegrationPoint{.weight = tet_qweights(n_points(), i),
                                .point = tet_qpoints(n_points(), i)};
    }
}
