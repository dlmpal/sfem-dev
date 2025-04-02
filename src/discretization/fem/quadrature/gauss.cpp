#include "gauss.hpp"

namespace sfem::fem::quadrature
{
    //=============================================================================
    static real_t gauss_1d_qweights(int n_points, int i)
    {
        static const real_t gauss_qwt_1[] = {2.0};

        static const real_t gauss_qwt_2[] = {1.0,
                                             1.0};

        static const real_t gauss_qwt_3[] = {5.0 / 9.0,
                                             8.0 / 9.0,
                                             5.0 / 9.0};

        static const real_t gauss_qwt_4[] = {0.347854845137454,
                                             0.652145154862546,
                                             0.652145154862546,
                                             0.347854845137454};

        static const real_t gauss_qwt_5[] = {0.236926885056189,
                                             0.478628670499366,
                                             0.568888888888889,
                                             0.478628670499366,
                                             0.236926885056189};

        switch (n_points)
        {
        case 1:
            return gauss_qwt_1[0];
        case 2:
            return gauss_qwt_2[i];
        case 3:
            return gauss_qwt_3[i];
        case 4:
            return gauss_qwt_4[i];
        case 5:
            return gauss_qwt_5[i];
        default:
            return 0;
        }
    }
    //=============================================================================
    static real_t gauss_1d_qpoints(int n_points, int i)
    {
        static const real_t gauss_qpt_1[] = {0.0};

        static const real_t gauss_qpt_2[] = {-0.577350269189626,
                                             0.577350269189626};

        static const real_t gauss_qpt_3[] = {-0.774596669241483,
                                             0.0,
                                             0.774596669241483};

        static const real_t gauss_qpt_4[] = {-0.861136311594053,
                                             -0.339981043584856,
                                             0.339981043584856,
                                             0.861136311594053};

        static const real_t gauss_qpt_5[] = {-0.906179845938664,
                                             -0.538469310105683,
                                             0.0,
                                             0.538469310105683,
                                             0.906179845938664};

        switch (n_points)
        {
        case 1:
            return gauss_qpt_1[0];
        case 2:
            return gauss_qpt_2[i];
        case 3:
            return gauss_qpt_3[i];
        case 4:
            return gauss_qpt_4[i];
        case 5:
            return gauss_qpt_5[i];
        default:
            return 0;
        }
    }
    //=============================================================================
    int Gauss1D::n_points() const
    {
        return order_ + 1;
    }
    //=============================================================================
    IntegrationPoint Gauss1D::point(int i) const
    {
        return IntegrationPoint{
            .weight = gauss_1d_qweights(order_ + 1, i),
            .point = {gauss_1d_qpoints(order_ + 1, i), 0, 0}};
    }
    //=============================================================================
    int Gauss2D::n_points() const
    {
        const int p = order_ + 1;
        return p * p;
    }
    //=============================================================================
    IntegrationPoint Gauss2D::point(int i) const
    {
        const int p = order_ + 1;
        return IntegrationPoint{
            .weight = gauss_1d_qweights(p, i % p) * gauss_1d_qweights(p, i / p),
            .point = {gauss_1d_qpoints(p, i % p), gauss_1d_qpoints(p, i / p)}};
    }
    //=============================================================================
    int Gauss3D::n_points() const
    {
        const int p = order_ + 1;
        return p * p * p;
    }
    //=============================================================================
    IntegrationPoint Gauss3D::point(int i) const
    {
        const int p = order_ + 1;
        return IntegrationPoint{
            .weight = gauss_1d_qweights(p, i % p) *
                      gauss_1d_qweights(p, (i % (p * p)) / p) *
                      gauss_1d_qweights(p, i / (p / p)),
            .point = {
                gauss_1d_qpoints(p, i % p),
                gauss_1d_qpoints(p, (i % (p * p)) / p),
                gauss_1d_qpoints(p, i / (p / p))}};
    }
}