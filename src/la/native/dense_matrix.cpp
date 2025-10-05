#include "dense_matrix.hpp"
#include "dense_matrix_utils.hpp"
#include "../../base/error.hpp"
#include <cmath>
#include <format>
#include <numeric>

namespace sfem::la
{
    //=============================================================================
    DenseMatrix::DenseMatrix(int n_rows, int n_cols, std::vector<real_t> &&values)
        : n_rows_(n_rows),
          n_cols_(n_cols),
          values_(std::move(values))
    {
        if (n_rows_ <= 0 || n_cols_ <= 0)
        {
            SFEM_ERROR(std::format("Invalid number of rows {}, or columns {}\n", n_rows_, n_cols_));
        }
        SFEM_CHECK_SIZES(n_rows_ * n_cols_, values_.size());
    }
    //=============================================================================
    DenseMatrix::DenseMatrix(int n_rows, int n_cols, real_t value)
        : DenseMatrix(n_rows, n_cols, std::vector<real_t>(n_rows * n_cols, value))
    {
    }
    //=============================================================================
    int DenseMatrix::n_rows() const
    {
        return n_rows_;
    }
    //=============================================================================
    int DenseMatrix::n_cols() const
    {
        return n_cols_;
    }
    //=============================================================================
    std::vector<real_t> &DenseMatrix::values()
    {
        return values_;
    }
    //=============================================================================
    const std::vector<real_t> &DenseMatrix::values() const
    {
        return values_;
    }
    //=============================================================================
    void DenseMatrix::set_all(real_t value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }
    //=============================================================================
    DenseMatrix DenseMatrix::copy() const
    {
        std::vector<real_t> data = values_;
        return DenseMatrix(n_rows_, n_cols_, std::move(data));
    }
    //=============================================================================
    DenseMatrix DenseMatrix::transpose() const
    {
        DenseMatrix transpose(n_cols_, n_rows_);
        utils::transpose(n_rows_, n_cols_, values_, transpose.values_);
        return transpose;
    }
    //=============================================================================
    DenseMatrix DenseMatrix::T() const
    {
        return transpose();
    }
    //=============================================================================
    std::pair<DenseMatrix, real_t> DenseMatrix::invert() const
    {
        DenseMatrix inv(n_cols_, n_rows_);
        real_t det;
        if (n_rows_ == n_cols_)
        {
            det = utils::inv(n_rows_, values_, inv.values_);
        }
        else
        {
            det = utils::pinv(n_rows_, n_cols_, values_, inv.values_);
        }
        return {std::move(inv), det};
    }
    //=============================================================================
    std::vector<real_t> DenseMatrix::row(int row_idx) const
    {
        std::vector<real_t> row(n_cols_);
        for (int i = 0; i < n_cols_; i++)
        {
            row[i] = values_[row_idx * n_cols_ + i];
        }
        return row;
    }
    //=============================================================================
    std::vector<real_t> DenseMatrix::col(int col_idx) const
    {
        std::vector<real_t> col(n_rows_);
        for (int i = 0; i < n_rows_; i++)
        {
            col[i] = values_[i * n_cols_ + col_idx];
        }
        return col;
    }
    //=============================================================================
    real_t &DenseMatrix::operator()(int i, int j)
    {
        /// @todo Add bounds check
        return values_[i * n_cols_ + j];
    }
    //=============================================================================
    real_t DenseMatrix::operator()(int i, int j) const
    {
        /// @todo Add bounds check
        return values_[i * n_cols_ + j];
    }
    //=============================================================================
    DenseMatrix &DenseMatrix::operator+=(real_t alpha)
    {
        for (std::size_t i = 0; i < values_.size(); i++)
        {
            values_[i] += alpha;
        }
        return *this;
    }
    //=============================================================================
    DenseMatrix &DenseMatrix::operator+=(const DenseMatrix &other)
    {
        SFEM_CHECK_SIZES(n_rows_, other.n_rows_);
        SFEM_CHECK_SIZES(n_cols_, other.n_cols_);
        utils::matadd(n_rows_, n_cols_, values_, 1, other.values_, 1, values_);
        return *this;
    }
    //=============================================================================
    DenseMatrix &DenseMatrix::operator-=(real_t alpha)
    {
        *this += -alpha;
        return *this;
    }
    //=============================================================================
    DenseMatrix &DenseMatrix::operator-=(const DenseMatrix &other)
    {
        SFEM_CHECK_SIZES(n_rows_, other.n_rows_);
        SFEM_CHECK_SIZES(n_cols_, other.n_cols_);
        utils::matadd(n_rows_, n_cols_, values_, 1, other.values_, -1, values_);
        return *this;
    }
    //=============================================================================
    DenseMatrix &DenseMatrix::operator*=(real_t alpha)
    {
        for (std::size_t i = 0; i < values_.size(); i++)
        {
            values_[i] *= alpha;
        }
        return *this;
    }
    //=============================================================================
    std::string DenseMatrix::str() const
    {
        std::string str = "[";
        for (int i = 0; i < n_rows_; i++)
        {
            if (i > 0)
            {
                str += " ";
            }
            for (int j = 0; j < n_cols_; j++)
            {
                str += std::format("{}", operator()(i, j));
                if (j < n_cols_ - 1)
                {
                    str += " ";
                }
            }
            if (i < n_rows_ - 1)
            {
                str += "\n";
            }
        }
        str += "]\n";
        return str;
    }
    //=============================================================================
    DenseMatrix operator+(const DenseMatrix &lhs, real_t rhs)
    {
        DenseMatrix result = lhs.copy();
        result += rhs;
        return result;
    }
    //=============================================================================
    DenseMatrix operator+(const DenseMatrix &lhs, const DenseMatrix &rhs)
    {
        DenseMatrix result = lhs.copy();
        result += rhs;
        return result;
    }
    //=============================================================================
    DenseMatrix operator-(const DenseMatrix &lhs, real_t rhs)
    {
        DenseMatrix result = lhs.copy();
        result -= rhs;
        return result;
    }
    //=============================================================================
    DenseMatrix operator-(const DenseMatrix &lhs, const DenseMatrix &rhs)
    {
        DenseMatrix result = lhs.copy();
        result += rhs;
        return result;
    }
    //=============================================================================
    DenseMatrix operator*(const DenseMatrix &lhs, const DenseMatrix &rhs)
    {
        SFEM_CHECK_SIZES(lhs.n_cols(), rhs.n_rows());
        DenseMatrix result(lhs.n_rows(), rhs.n_cols());
        utils::matmult(result.n_rows(), result.n_cols(), lhs.n_cols(),
                       lhs.values(), rhs.values(), result.values());
        return result;
    }
    //=============================================================================
    DenseMatrix operator*(const DenseMatrix &lhs, real_t rhs)
    {
        DenseMatrix result = lhs.copy();
        result *= rhs;
        return result;
    }
    //=============================================================================
    DenseMatrix submatrix(const DenseMatrix &mat,
                          int start_row, int end_row,
                          int start_col, int end_col)
    {
        DenseMatrix sub(end_row - start_row,
                        end_col - start_col);
        for (int i = start_row; i < end_row; i++)
        {
            for (int j = start_col; j < end_col; j++)
            {
                sub(i - start_row, j - start_col) = mat(i, j);
            }
        }
        return sub;
    }
    //=============================================================================
    real_t norm(const DenseMatrix &A)
    {
        return std::sqrt(std::accumulate(A.values().cbegin(),
                                         A.values().cend(), 0.0,
                                         [](real_t acc, real_t v)
                                         {
                                             return acc + v * v;
                                         }));
    }
}