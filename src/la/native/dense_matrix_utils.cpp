#include "dense_matrix_utils.hpp"
#include "../../base/error.hpp"
#include <cmath>
#include <vector>
#include <limits>

namespace sfem::la::utils
{
    //=============================================================================
    void transpose(int nr, int nc, std::span<const real_t> m, std::span<real_t> mt)
    {
        for (int i = 0; i < nr; i++)
        {
            for (int j = 0; j < nc; j++)
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
    real_t det3x3(int nr, std::span<const real_t> m)
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
    real_t inv3x3(int nr, std::span<const real_t> m, std::span<real_t> mi)
    {
        // Determinant and inverse
        real_t d = det3x3(nr, m);
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
    real_t inv(int nr, std::span<const real_t> A, std::span<real_t> Ainv)
    {
        if (nr <= 3)
        {
            return inv3x3(nr, A, Ainv);
        }

        // Matrix copy
        std::vector<real_t> A_(A.data(), A.data() + nr * nr);

        // Initialize inverse as identity
        for (int i = 0; i < nr; i++)
        {
            Ainv[i * nr + i] = 1.0;
        }

        const real_t eps = std::numeric_limits<real_t>::epsilon();

        // Swap the rows r1 and r2
        auto row_swap = [nr](std::span<real_t> mat, int r1, int r2)
        {
            for (int j = 0; j < nr; j++)
            {
                std::swap(mat[r1 * nr + j], mat[r2 * nr + j]);
            }
        };

        // Scale row r1 by a
        auto row_scale = [nr](std::span<real_t> mat, int r, real_t a)
        {
            for (int j = 0; j < nr; j++)
            {
                mat[r * nr + j] *= a;
            }
        };

        // Add row r1 (scaled by a) to row r2
        auto row_axpy = [nr](std::span<real_t> mat, int r1, int r2, real_t a)
        {
            for (int j = 0; j < nr; j++)
            {
                mat[r2 * nr + j] += a * mat[r1 * nr + j];
            }
        };

        int n_swap = 0;
        real_t det = 1.0;
        for (int c = 0; c < nr; c++)
        {
            // Find pivot row
            int p = c;
            for (int r = c; r < nr; r++)
            {
                if (std::abs(A_[r * nr + c]) > std::abs(A_[p * nr + c]))
                {
                    p = r;
                }
            }

            // Store pivot and check for singularity
            const real_t pivot = A_[p * nr + c];
            if (std::abs(pivot) < eps)
            {
                SFEM_ERROR("Singular matrix!\n");
            }

            det *= pivot;

            // Swap rows
            if (c != p)
            {
                row_swap(A_, c, p);
                row_swap(Ainv, c, p);
                n_swap++;
            }

            // Normalize pivot row
            row_scale(A_, c, 1 / pivot);
            row_scale(Ainv, c, 1 / pivot);

            // Elimininate rows
            for (int r = 0; r < nr; r++)
            {
                if (r == c)
                {
                    continue;
                }

                const real_t factor = -A_[r * nr + c];
                row_axpy(A_, c, r, factor);
                row_axpy(Ainv, c, r, factor);
            }
        }

        // Adjust determinant's sign based on number of row swaps
        if ((n_swap - 1) % 2 == 0)
        {
            det = -det;
        }

        return det;
    }
    //=============================================================================
    real_t pinv3x3(int nr, int nc, std::span<const real_t> m, std::span<real_t> mi)
    {
        // Transpose
        std::array<real_t, 3 * 3> mt;
        transpose(nr, nc, m, mt);

        // Matrix-transpose product
        std::array<real_t, 3 * 3> mtm;
        matmult(nc, nc, nr, mt, m, mtm);

        // Invert matrix-transpose product
        std::array<real_t, 3 * 3> mtmi;
        const real_t det = std::sqrt(inv3x3(nc, mtm, mtmi));

        // Multiply by tranpose
        matmult(nc, nr, nc, mtmi, mt, mi);

        return det;
    }
    //=============================================================================
    real_t pinv(int nr, int nc, std::span<const real_t> m, std::span<real_t> mi)
    {
        if (nr < nc)
        {
            SFEM_ERROR(std::format("Cannot compute pseudo-inverse for matrix with more columns ({}) than rows ({})\n", nc, nr));
        }

        if (nr <= 3)
        {
            return pinv3x3(nr, nc, m, mi);
        }

        // Transpose
        std::vector<real_t> mt(nc * nr);
        transpose(nr, nc, m, mt);

        // Matrix-transpose product
        std::vector<real_t> mtm(nc * nc);
        matmult(nc, nc, nr, mt, m, mtm);

        // Invert matrix-transpose product
        std::vector<real_t> mtmi(nc * nc);
        const real_t det = std::sqrt(inv(nc, mtm, mtmi));

        // Multiply by tranpose
        matmult(nc, nr, nc, mtmi, mt, mi);

        return det;
    }
}