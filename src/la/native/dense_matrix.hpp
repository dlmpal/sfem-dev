#pragma once

#include "../../base/config.hpp"
#include <array>
#include <vector>
#include <span>
#include <string>

namespace sfem::la
{
    class DenseMatrix
    {
    public:
        DenseMatrix(int n_rows, int n_cols, std::vector<real_t> &&values);

        DenseMatrix(int n_rows, int n_cols, real_t value = 0.0);

        DenseMatrix(const DenseMatrix &) = delete;
        DenseMatrix &operator=(const DenseMatrix &) = delete;

        DenseMatrix(DenseMatrix &&) = default;
        DenseMatrix &operator=(DenseMatrix &&) = default;

        int n_rows() const;

        int n_cols() const;

        std::vector<real_t> &data();

        const std::vector<real_t> &data() const;

        DenseMatrix copy() const;

        DenseMatrix transpose() const;

        DenseMatrix T() const;

        std::pair<DenseMatrix, real_t> invert() const;

        DenseMatrix row(int row_idx) const;

        DenseMatrix col(int col_idx) const;

        real_t operator()(int i, int j) const;

        real_t &operator()(int i, int j);

        DenseMatrix &operator+=(real_t alpha);

        DenseMatrix &operator+=(const DenseMatrix &other);

        DenseMatrix &operator-=(real_t alpha);

        DenseMatrix &operator-=(const DenseMatrix &other);

        DenseMatrix &operator*=(real_t alpha);

        std::string str() const;

    private:
        int n_rows_;

        int n_cols_;

        std::vector<real_t> data_;
    };

    DenseMatrix operator+(const DenseMatrix &lhs, real_t rhs);

    DenseMatrix operator+(const DenseMatrix &lhs, const DenseMatrix &rhs);

    DenseMatrix operator-(const DenseMatrix &lhs, real_t rhs);

    DenseMatrix operator-(const DenseMatrix &lhs, const DenseMatrix &rhs);

    DenseMatrix operator*(const DenseMatrix &lhs, const DenseMatrix &rhs);

    DenseMatrix operator*(const DenseMatrix &lhs, real_t rhs);
}