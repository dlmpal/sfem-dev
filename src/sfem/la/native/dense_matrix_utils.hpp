#pragma once

#include <sfem/base/config.hpp>
#include <span>

namespace sfem::la::utils
{
    /// @brief Tranpose a matrix
    /// @param nr Number of rows
    /// @param nc Number of columns
    /// @param m Matrix
    /// @param mt Matrix transpose
    void transpose(int nr, int nc,
                   std::span<const real_t> m,
                   std::span<real_t> mt);

    /// @brief Matrix multiplication: m3 = m1 * m2
    /// @param nr Number of rows of first matrix
    /// @param nc Number of columns of second matrix
    /// @param nrc Number of columns of first matrix (=number of columns of second)
    /// @param m1 First matrix
    /// @param m2 Second matrix
    /// @param m3 Result
    void matmult(int nr, int nc, int nrc,
                 std::span<const real_t> m1,
                 std::span<const real_t> m2,
                 std::span<real_t> m3);

    /// @brief Scaled matrix addition: m3 = a1 * m1 + a2 * m2
    /// @param nr Number of rows
    /// @param nc Number of columns
    /// @param m1 First matrix
    /// @param a1 First scalar
    /// @param m2 Second matrix
    /// @param a2 Second scalar
    /// @param m3 Result
    void matadd(int nr, int nc,
                std::span<const real_t> m1, real_t a1,
                std::span<const real_t> m2, real_t a2,
                std::span<real_t> m3);

    /// @brief Compute the inverse and determinant of square matrix
    /// @param nr Number of rows (=number of columns)
    /// @param m Matrix
    /// @param mi Inverse matrix
    /// @return Determinant
    /// @note Optimized for 3x3, 2x2 and 1x1 matrices
    real_t inv(int r, std::span<const real_t> m, std::span<real_t> mi);

    /// @brief Compute the Moore-Penrose pseudo-inverse of a matrix (and its determinant)
    /// @param nr Number of rows
    /// @param nc Number of columns
    /// @param m Matrix
    /// @param mi Matrix pseudo-inverse
    /// @return Determinant
    /// @note Optimized for 3x2, 3x1 and 2x1 matrices
    real_t pinv(int r, int c, std::span<const real_t> m, std::span<real_t> mi);
}
