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

        static const real_t tet_qwt_14[] = {0.0070910034628469025, 0.007091003462846909, 0.007091003462846909,
                                            0.007091003462846912, 0.007091003462846912, 0.0070910034628469155,
                                            0.012248840519393652, 0.012248840519393652, 0.012248840519393655,
                                            0.012248840519393659, 0.018781320953002632, 0.018781320953002632,
                                            0.018781320953002632, 0.01878132095300265};

        switch (n_points)
        {
        case 4:
            return tet_qwt_4[i];
        case 5:
            return tet_qwt_5[i];
        case 14:
            return tet_qwt_14[i];
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

        static const real_t tet_qpt_14[] = {0.4544962958743503, 0.4544962958743504, 0.04550370412564962,
                                            0.04550370412564967, 0.4544962958743504, 0.4544962958743504,
                                            0.04550370412564973, 0.4544962958743503, 0.04550370412564969,
                                            0.4544962958743503, 0.04550370412564966, 0.4544962958743504,
                                            0.4544962958743503, 0.04550370412564968, 0.04550370412564962,
                                            0.0455037041256497, 0.04550370412564966, 0.4544962958743504,
                                            0.09273525031089128, 0.7217942490673263, 0.09273525031089122,
                                            0.721794249067326, 0.09273525031089128, 0.09273525031089129,
                                            0.09273525031089132, 0.09273525031089114, 0.09273525031089129,
                                            0.0927352503108913, 0.0927352503108913, 0.7217942490673263,
                                            0.3108859192633006, 0.06734224221009831, 0.3108859192633006,
                                            0.06734224221009824, 0.3108859192633006, 0.3108859192633007,
                                            0.3108859192633006, 0.3108859192633007, 0.3108859192633006,
                                            0.3108859192633006, 0.3108859192633007, 0.06734224221009814};
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
        case 14:
            return {tet_qpt_14[i * 3 + 0],
                    tet_qpt_14[i * 3 + 1],
                    tet_qpt_14[i * 3 + 2]};
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
