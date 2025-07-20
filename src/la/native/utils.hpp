#pragma once

#include "../../base/error.hpp"
#include "../../base/config.hpp"
#include <cmath>
#include <vector>

namespace sfem::la::utils
{
    /// @brief Tranpose a matrix of size (r, c)
    inline void transpose(int r, int c,
                          const real_t m[], real_t mt[])
    {
        for (auto i = 0; i < r; i++)
        {
            for (auto j = 0; j < c; j++)
            {
                mt[j * r + i] = m[i * c + j];
            }
        }
    }

    /// @brief Multiply two matrices of size (r, rc) and (rc, c)
    inline void matmult(int r, int c, int rc,
                        const real_t m1[], const real_t m2[],
                        real_t m[])
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                m[i * c + j] = 0.0;

                for (int k = 0; k < rc; k++)
                {
                    m[i * c + j] = std::fma(m1[i * rc + k], m2[k * c + j], m[i * c + j]);
                }
            }
        }
    }

    /// @brief Add two matrices of size (r, c)
    /// @note the matrix are scaled by a1 and a2, i.e. m = a1 * m1 + a2 * m2
    inline void matadd(int r, int c,
                       const real_t m1[], real_t a1,
                       const real_t m2[], real_t a2,
                       real_t m[])
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                m[i * c + j] = a1 * m1[i * c + j] + a2 * m2[i * c + j];
            }
        }
    }

    /// @brief Compute y += a*x
    inline void axpy(int size, real_t a, const real_t x[], real_t y[])
    {
        for (int i = 0; i < size; i++)
        {
            y[i] = std::fma(a, x[i], y[i]);
        }
    }

    /// @brief Determinant of a square 3x3, 2x2 or 1x1 matrix
    inline real_t det(int r, const real_t m[])
    {
        real_t d = -1;

        if (r == 3)
        {
            d = m[0] * (m[4] * m[8] - m[7] * m[5]) -
                m[1] * (m[3] * m[8] - m[6] * m[5]) +
                m[2] * (m[3] * m[7] - m[4] * m[6]);
        }
        else if (r == 2)
        {
            d = m[0] * m[3] - m[1] * m[2];
        }
        else if (r == 1)
        {
            d = m[0];
        }
        else
        {
            SFEM_ERROR("Determinant not defined for square matrix with dimension: " + std::to_string(r) + "\n");
        }

        return d;
    }

    /// @brief Inverse and determinant of a square 3x3, 2x2 or 1x1 matrix
    inline real_t inv(int r, const real_t m[], real_t mi[])
    {
        // Determinant and inverse
        real_t d = det(r, m);
        real_t di = 1 / d;

        if (r == 3)
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
        else if (r == 2)
        {
            mi[0] = di * m[3];
            mi[1] = -di * m[1];
            mi[2] = -di * m[2];
            mi[3] = di * m[0];
        }
        else if (r == 1)
        {
            mi[0] = di;
        }
        else
        {
            SFEM_ERROR("Inverse not defined for square matrix with dimension: " + std::to_string(r) + "\n");
        }

        return d;
    }

    /// @brief Moore-Penrose pseudo-inverse of a 3x2, 3x1 or 2x1 matrix
    inline real_t pinv(int r, int c, const real_t m[], real_t mi[])
    {
        // Transpose
        real_t mt[3 * 3];
        transpose(r, c, m, mt);

        // Intermediate product
        real_t mtm[3 * 3];
        matmult(c, c, r, mt, m, mtm);

        // Invert intermediate product
        real_t mtmi[3 * 3] = {}; /// < Suppress unitialized warning
        inv(c, mtm, mtmi);

        // Multiply by tranpose
        matmult(c, r, c, mtmi, mt, mi);

        // Compute the determinant
        real_t d = sqrt(det(c, mtm));

        return d;
    }
}
