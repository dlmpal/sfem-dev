#pragma once

#include "../../base/config.hpp"
#include <array>
#include <vector>
#include <span>
#include <string>

namespace sfem::la
{
    /// @brief Dense matrix
    class DenseMatrix
    {
    public:
        /// @brief Create a DenseMatrix
        /// @param n_rows Number of rows
        /// @param n_cols Number of columns
        /// @param values Matrix values
        DenseMatrix(int n_rows, int n_cols, std::vector<real_t> &&values);

        /// @brief Create a DenseMatrix
        /// @param n_rows Number of rows
        /// @param n_cols Number of columns
        /// @param values Uniform value
        DenseMatrix(int n_rows, int n_cols, real_t value = 0.0);

        // Avoid unintentional copies
        DenseMatrix(const DenseMatrix &) = delete;
        DenseMatrix &operator=(const DenseMatrix &) = delete;

        // Move constructor and assignment
        DenseMatrix(DenseMatrix &&) = default;
        DenseMatrix &operator=(DenseMatrix &&) = default;

        /// @brief Get the number of rows
        int n_rows() const;

        /// @brief Get the number of columns
        int n_cols() const;

        /// @brief Get the matrix values
        std::vector<real_t> &values();

        /// @brief Get the matrix values (const version)
        const std::vector<real_t> &values() const;

        /// @brief Set all matrix values to a uniform value
        void set_all(real_t value);

        /// @brief Get a copy of the matrix
        DenseMatrix copy() const;

        /// @brief Get the transpose of the matrix
        DenseMatrix transpose() const;

        /// @brief Get the transpose of the matrix
        DenseMatrix T() const;

        /// @brief Invert the matrix
        /// @note Works only for 1x1, 2x2 and 3x3 matrices
        std::pair<DenseMatrix, real_t> invert() const;

        /// @brief Get the values of a given row
        std::vector<real_t> row(int row_idx) const;

        /// @brief Get the values of a given column
        std::vector<real_t> col(int col_idx) const;

        /// @brief Get the value at a given index pair
        real_t &operator()(int i, int j);

        /// @brief Get the value at a given index pair (const version)
        real_t operator()(int i, int j) const;

        // Arithmetic operations
        DenseMatrix &operator+=(real_t alpha);
        DenseMatrix &operator+=(const DenseMatrix &other);
        DenseMatrix &operator-=(real_t alpha);
        DenseMatrix &operator-=(const DenseMatrix &other);
        DenseMatrix &operator*=(real_t alpha);

        /// @brief Get a string representation of the matrix
        std::string str() const;

    private:
        /// @brief Number of rows
        int n_rows_;

        /// @brief Number of columns
        int n_cols_;

        /// @brief Values
        std::vector<real_t> values_;
    };

    // Arithmetic operators
    DenseMatrix operator+(const DenseMatrix &lhs, real_t rhs);
    DenseMatrix operator+(const DenseMatrix &lhs, const DenseMatrix &rhs);
    DenseMatrix operator-(const DenseMatrix &lhs, real_t rhs);
    DenseMatrix operator-(const DenseMatrix &lhs, const DenseMatrix &rhs);
    DenseMatrix operator*(const DenseMatrix &lhs, const DenseMatrix &rhs);
    DenseMatrix operator*(const DenseMatrix &lhs, real_t rhs);

    /// @brief Extract a submatrix
    DenseMatrix submatrix(const DenseMatrix &mat,
                          int start_row, int end_row,
                          int start_col, int end_col);

    /// @brief Compute the Frobenius norm for a dense matrix
    real_t norm(const DenseMatrix &A);
}