#include "dense_matrix_utils.hpp"
#include "../../base/error.hpp"
#include <cmath>

namespace sfem::la::utils
{
    //=============================================================================
    void transpose(int nr, int nc, std::span<const real_t> m, std::span<real_t> mt)
    {
        for (auto i = 0; i < nr; i++)
        {
            for (auto j = 0; j < nc; j++)
            {
                mt[j * nr + i] = m[i * nc + j];
            }
        }
    }
    //=============================================================================
    void matmult(int nr, int nc, int nrc,
                 std::span<const real_t> m1,
                 std::span<const real_t> m2,
                 std::span<real_t> m3)
    {
        for (int i = 0; i < nr; i++)
        {
            for (int j = 0; j < nc; j++)
            {
                m3[i * nc + j] = 0.0;

                for (int k = 0; k < nrc; k++)
                {
                    m3[i * nc + j] = std::fma(m1[i * nrc + k],
                                              m2[k * nc + j],
                                              m3[i * nc + j]);
                }
            }
        }
    }
    //=============================================================================
    void matadd(int nr, int nc,
                std::span<const real_t> m1, real_t a1,
                std::span<const real_t> m2, real_t a2,
                std::span<real_t> m3)
    {
        for (int i = 0; i < nr; i++)
        {
            for (int j = 0; j < nc; j++)
            {
                m3[i * nc + j] = a1 * m1[i * nc + j] + a2 * m2[i * nc + j];
            }
        }
    }
    //=============================================================================
    real_t det(int nr, std::span<const real_t> m)
    {
        real_t d = -1;

        if (nr == 3)
        {
            d = m[0] * (m[4] * m[8] - m[7] * m[5]) -
                m[1] * (m[3] * m[8] - m[6] * m[5]) +
                m[2] * (m[3] * m[7] - m[4] * m[6]);
        }
        else if (nr == 2)
        {
            d = m[0] * m[3] - m[1] * m[2];
        }
        else if (nr == 1)
        {
            d = m[0];
        }
        else
        {
            SFEM_ERROR("Determinant not defined for square matrix with dimension: " + std::to_string(nr) + "\n");
        }

        return d;
    }
    //=============================================================================
    real_t inv(int nr, std::span<const real_t> m, std::span<real_t> mi)
    {
        // Determinant and inverse
        real_t d = det(nr, m);
        real_t di = 1 / d;

        if (nr == 3)
        {
            mi[0] = di * (m[4] * m[8] - m[5] * m[7]);
            mi[1] = -di * (m[1] * m[8] - m[2] * m[7]);
            mi[2] = di * (m[1] * m[5] - m[2] * m[4]);
            mi[3] = -di * (m[3] * m[8] - m[5] * m[6]);
            mi[4] = di * (m[0] * m[8] - m[2] * m[6]);
            mi[5] = -di * (m[0] * m[5] - m[2] * m[3]);
            mi[6] = di * (m[3] * m[7] - m[4] * m[6]);
            mi[7] = -di * (m[0] * m[7] - m[1] * m[6]);
            mi[8] = di * (m[0] * m[4] - m[1] * m[3]);
        }
        else if (nr == 2)
        {
            mi[0] = di * m[3];
            mi[1] = -di * m[1];
            mi[2] = -di * m[2];
            mi[3] = di * m[0];
        }
        else if (nr == 1)
        {
            mi[0] = di;
        }
        else
        {
            SFEM_ERROR("Inverse not defined for square matrix with dimension: " + std::to_string(nr) + "\n");
        }

        return d;
    }
    //=============================================================================
    real_t pinv(int nr, int nc, std::span<const real_t> m, std::span<real_t> mi)
    {
        // Transpose
        std::array<real_t, 3 * 3> mt;
        transpose(nr, nc, m, mt);

        // Intermediate product
        std::array<real_t, 3 * 3> mtm;
        matmult(nc, nc, nr, mt, m, mtm);

        // Invert intermediate product
        std::array<real_t, 3 * 3> mtmi; /// < Suppress unitialized warning
        inv(nc, mtm, mtmi);

        // Multiply by tranpose
        matmult(nc, nr, nc, mtmi, mt, mi);

        // Compute the determinant
        real_t d = std::sqrt(det(nc, mtm));

        return d;
    }
}