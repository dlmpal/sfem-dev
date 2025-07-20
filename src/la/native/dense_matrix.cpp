#include "dense_matrix.hpp"
#include "utils.hpp"
#include <format>

namespace sfem::la
{
    //=============================================================================
    DenseMatrix::DenseMatrix(int n_rows, int n_cols, std::vector<real_t> &&values)
        : n_rows_(n_rows),
          n_cols_(n_cols),
          values_(values)
    {
        // Check that sizes match
        if (n_rows_ <= 0 || n_cols_ <= 0)
        {
            /// @todo error
        }
        if (n_rows_ * n_cols_ != static_cast<int>(values.size()))
        {
            /// @todo error
        }
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
    DenseMatrix DenseMatrix::copy() const
    {
        std::vector<real_t> data = values_;
        return DenseMatrix(n_rows_, n_cols_, std::move(data));
    }
    //=============================================================================
    DenseMatrix DenseMatrix::transpose() const
    {
        DenseMatrix transpose(n_cols_, n_rows_);
        utils::transpose(n_rows_,
                         n_cols_,
                         values_.data(),
                         transpose.values_.data());
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
            det = utils::inv(n_rows_,
                             values_.data(),
                             inv.values_.data());
        }
        else
        {
            det = utils::pinv(n_rows_,
                              n_cols_,
                              values_.data(),
                              inv.values_.data());
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
        return values_[i * n_cols_ + j];
    }
    //=============================================================================
    real_t DenseMatrix::operator()(int i, int j) const
    {
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
        utils::matadd(n_rows_, n_cols_,
                      values_.data(), 1,
                      other.values_.data(), 1,
                      values_.data());
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
        utils::matadd(n_rows_, n_cols_,
                      values_.data(), 1,
                      other.values_.data(), -1,
                      values_.data());
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
        DenseMatrix result(lhs.n_rows(), rhs.n_cols());
        utils::matadd(result.n_rows(),
                      result.n_cols(),
                      lhs.values().data(), 1,
                      rhs.values().data(), 1,
                      result.values().data());
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
        DenseMatrix result(lhs.n_rows(), rhs.n_cols());
        utils::matadd(result.n_rows(),
                      result.n_cols(),
                      lhs.values().data(), 1,
                      rhs.values().data(), -1,
                      result.values().data());
        return result;
    }
    //=============================================================================
    DenseMatrix operator*(const DenseMatrix &lhs, const DenseMatrix &rhs)
    {
        DenseMatrix result(lhs.n_rows(), rhs.n_cols());
        utils::matmult(result.n_rows(),
                       result.n_cols(),
                       lhs.n_cols(),
                       lhs.values().data(),
                       rhs.values().data(),
                       result.values().data());
        return result;
    }
    //=============================================================================
    DenseMatrix operator*(const DenseMatrix &lhs, real_t rhs)
    {
        DenseMatrix result = lhs.copy();
        result *= rhs;
        return result;
    }
}