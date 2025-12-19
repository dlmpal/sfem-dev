#include "vector.hpp"
#include <sfem/parallel/mpi.hpp>
#include <sfem/parallel/scatterer.hpp>
#include <sfem/base/error.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace sfem::la
{
    //=============================================================================
    Vector::Vector(std::shared_ptr<const IndexMap> im,
                   int block_size,
                   std::vector<real_t> &&values)
        : im_(im),
          bs_(block_size),
          values_(std::move(values))
    {
        SFEM_CHECK_SIZES(im_->n_local() * bs_, values_.size());
    }
    //=============================================================================
    Vector::Vector(std::shared_ptr<const IndexMap> im,
                   int block_size,
                   real_t value)
        : Vector(im, block_size,
                 std::vector<real_t>(im->n_local() * block_size, value))
    {
    }
    //=============================================================================
    std::shared_ptr<const IndexMap> Vector::index_map() const
    {
        return im_;
    }
    //=============================================================================
    int Vector::block_size() const
    {
        return bs_;
    }
    //=============================================================================
    std::vector<real_t> &Vector::values()
    {
        return values_;
    }
    //=============================================================================
    const std::vector<real_t> &Vector::values() const
    {
        return values_;
    }
    //=============================================================================
    int Vector::n_owned() const
    {
        return im_->n_owned();
    }
    //=============================================================================
    int Vector::n_ghost() const
    {
        return im_->n_ghost();
    }
    //=============================================================================
    int Vector::n_local() const
    {
        return im_->n_local();
    }
    //=============================================================================
    int Vector::n_global() const
    {
        return im_->n_global();
    }
    //=============================================================================
    real_t &Vector::operator()(int idx, int comp)
    {
        return values_[idx * bs_ + comp];
    }
    //=============================================================================
    real_t Vector::operator()(int idx, int comp) const
    {
        return values_[idx * bs_ + comp];
    }
    //=============================================================================
    void Vector::set_all(real_t value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }
    //=============================================================================
    void Vector::set_values(std::span<const int> idxs,
                            std::span<const real_t> values,
                            SetMode mode)
    {
        SFEM_CHECK_SIZES(idxs.size() * bs_, values.size());
        if (mode == SetMode::insert)
        {
            for (std::size_t i = 0; i < idxs.size(); i++)
            {
                for (int j = 0; j < bs_; j++)
                {
                    values_[idxs[i] * bs_ + j] = values[i * bs_ + j];
                }
            }
        }
        else // SetMode::add
        {
            for (std::size_t i = 0; i < idxs.size(); i++)
            {
                for (int j = 0; j < bs_; j++)
                {
                    values_[idxs[i] * bs_ + j] += values[i * bs_ + j];
                }
            }
        }
    }
    //=============================================================================
    void Vector::assemble()
    {
        /// @todo Construct this only once
        Scatterer<real_t> scatter(im_);
        scatter.reverse(values_, bs_,
                        [](real_t &dest, real_t src)
                        { dest += src; });
        // Set ghost values to zero
        std::fill(values_.begin() + n_owned() * bs_, values_.end(), 0.0);
    }
    //=============================================================================
    void Vector::update_ghosts()
    {
        /// @todo Construct this only once
        Scatterer<real_t> scatter(im_);
        scatter.forward(values_, bs_,
                        [](real_t &dest, real_t src)
                        { dest = src; });
    }
    //=============================================================================
    void copy(const Vector &src, Vector &dest)
    {
        SFEM_CHECK_SIZES(src.block_size(), dest.block_size());
        SFEM_CHECK_SIZES(src.n_local(), dest.n_local());
        std::copy(src.values().cbegin(),
                  src.values().cbegin() + src.n_owned() * src.block_size(),
                  dest.values().begin());
    }
    //=============================================================================
    void scale(real_t a, Vector &x)
    {
        for (int i = 0; i < x.n_owned(); i++)
        {
            for (int j = 0; j < x.block_size(); j++)
            {
                x(i, j) *= a;
            }
        }
    }
    //=============================================================================
    void axpy(real_t a, const Vector &x, Vector &y)
    {
        SFEM_CHECK_SIZES(x.block_size(), y.block_size());
        SFEM_CHECK_SIZES(x.n_owned(), y.n_owned());
        for (int i = 0; i < x.n_owned(); i++)
        {
            for (int j = 0; j < x.block_size(); j++)
            {
                y(i, j) += a * x(i, j);
            }
        }
    }
    //=============================================================================
    void axpbypc(real_t a, real_t b, real_t c,
                 const Vector &x, const Vector &y, Vector &z)
    {
        SFEM_CHECK_SIZES(x.n_owned(), y.n_owned());
        SFEM_CHECK_SIZES(x.n_owned(), z.n_owned());
        SFEM_CHECK_SIZES(x.block_size(), y.block_size());
        SFEM_CHECK_SIZES(x.block_size(), z.block_size());
        for (int i = 0; i < x.n_owned(); i++)
        {
            for (int j = 0; j < x.block_size(); j++)
            {
                real_t val = a * x(i, j) + b * y(i, j) + c;
                z(i, j) = val;
            }
        }
    }
    //=============================================================================
    real_t dot(const Vector &x, const Vector &y)
    {
        SFEM_CHECK_SIZES(x.block_size(), y.block_size());
        SFEM_CHECK_SIZES(x.n_owned(), y.n_owned());
        real_t prod = std::inner_product(x.values().cbegin(),
                                         x.values().cbegin() + x.n_owned() * x.block_size(),
                                         y.values().cbegin(), 0.0);
        return mpi::reduce(prod, mpi::ReduceOperation::sum);
    }
    //=============================================================================
    real_t norm(const Vector &x, NormType norm_type)
    {
        real_t val = 0;
        switch (norm_type)
        {
        case NormType::l1:
            val = std::accumulate(x.values().cbegin(),
                                  x.values().cend(), 0.0,
                                  [](real_t acc, real_t v)
                                  {
                                      return acc + std::abs(v);
                                  });
            val = mpi::reduce(val, mpi::ReduceOperation::sum);
            break;
        case NormType::l2:
            val = std::sqrt(dot(x, x));
            break;
        case NormType::linf:
            val = *std::max_element(x.values().cbegin(),
                                    x.values().cend());
            val = mpi::reduce(val, mpi::ReduceOperation::max);
            break;
        default:
            SFEM_ERROR("Invalid norm type\n");
            break;
        }
        return val;
    }
}