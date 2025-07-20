#include "vector.hpp"
#include "../../parallel/mpi.hpp"
#include "../../base/error.hpp"
#include <fstream>
namespace sfem::la
{
    //=============================================================================
    Vector::Vector(std::shared_ptr<const IndexMap> im,
                   int block_size,
                   std::vector<real_t> &&values)
        : index_map_(im),
          block_size_(block_size),
          values_(std::move(values))
    {
        SFEM_CHECK_SIZES(index_map_->n_local() * block_size_, values_.size());
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
        return index_map_;
    }
    //=============================================================================
    int Vector::block_size() const
    {
        return block_size_;
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
    int Vector::n_local() const
    {
        return index_map_->n_local();
    }
    //=============================================================================
    int Vector::n_global() const
    {
        return index_map_->n_global();
    }
    //=============================================================================
    real_t &Vector::operator()(int idx, int comp)
    {
        return values_[idx * block_size_ + comp];
    }
    //=============================================================================
    real_t Vector::operator()(int idx, int comp) const
    {
        return values_[idx * block_size_ + comp];
    }
    //=============================================================================
    void Vector::set_all(real_t value)
    {
        std::fill(values_.begin(),
                  values_.begin() + index_map_->n_owned() * block_size_,
                  value);
        assemble();
    }
    //=============================================================================
    void Vector::set_values(std::span<const int> idxs,
                            std::span<const real_t> values,
                            bool insert)
    {
        SFEM_CHECK_SIZES(idxs.size() * block_size_, values.size());
        int bs = block_size_;
        if (insert)
        {
            for (std::size_t i = 0; i < idxs.size(); i++)
            {
                for (int j = 0; j < bs; j++)
                {
                    values_[idxs[i] * bs + j] = values[i * bs + j];
                }
            }
        }
        else
        {
            for (std::size_t i = 0; i < idxs.size(); i++)
            {
                for (int j = 0; j < bs; j++)
                {
                    values_[idxs[i] * bs + j] += values[i * bs + j];
                }
            }
        }
    }
    //=============================================================================
    void Vector::assemble()
    {
        int bs = block_size_;

        // Send ghost value contributions to owner processes
        auto [recv_ghost_idxs,
              recv_ghost_counts,
              recv_ghost_displs] = mpi::send_to_dest<int>(index_map_->ghost_idxs(),
                                                          index_map_->ghost_owners());
        auto recv_ghost_values = std::get<0>(mpi::send_to_dest<real_t>({values_.cbegin() + index_map_->n_owned() * bs,
                                                                        values_.cend()},
                                                                       index_map_->ghost_owners(), bs));

        // Add ghost value contributions
        recv_ghost_idxs = index_map_->global_to_local(recv_ghost_idxs);
        for (std::size_t i = 0; i < recv_ghost_idxs.size(); i++)
        {
            for (int j = 0; j < bs; j++)
            {
                values_[recv_ghost_idxs[i] * bs + j] += recv_ghost_values[i * bs + j];
            }
        }

        // After adding all the contributions, each process has to send
        // the recently computed values for indices for which it received
        // contributions back. This way, the ghost values for each process
        // are correct, and in agreement with the value in the owner process
        std::vector<int> send_ghost_dest(recv_ghost_idxs.size());
        std::vector<int> send_ghost_idxs(recv_ghost_idxs.size());
        std::vector<real_t> send_ghost_values(recv_ghost_idxs.size() * bs);
        for (int i = 0; i < mpi::n_procs(); i++)
        {
            for (int j = 0; j < recv_ghost_counts[i]; j++)
            {
                int pos = recv_ghost_displs[i] + j;
                send_ghost_dest[pos] = i;
                send_ghost_idxs[pos] = index_map_->local_to_global(recv_ghost_idxs[pos]);
                for (int k = 0; k < bs; k++)
                {
                    send_ghost_values[pos * bs + k] = values_[recv_ghost_idxs[pos] * bs + k];
                }
            }
        }
        std::tie(recv_ghost_idxs,
                 recv_ghost_counts,
                 recv_ghost_displs) = mpi::send_to_dest<int>(send_ghost_idxs, send_ghost_dest);
        recv_ghost_values = std::get<0>(mpi::send_to_dest<real_t>(send_ghost_values, send_ghost_dest, bs));

        // Update ghost values
        recv_ghost_idxs = index_map_->global_to_local(recv_ghost_idxs);
        for (int i = 0; i < index_map_->n_ghost(); i++)
        {
            for (int j = 0; j < bs; j++)
            {
                values_[recv_ghost_idxs[i] * bs + j] = recv_ghost_values[i * bs + j];
            }
        }
    }
    //=============================================================================
    void vec_copy(const Vector &src, Vector &dest)
    {
        SFEM_CHECK_SIZES(src.block_size(), dest.block_size());
        SFEM_CHECK_SIZES(src.n_local(), dest.n_local());
        const auto &src_values = src.values();
        auto &dest_values = dest.values();
        std::copy(src_values.cbegin(),
                  src_values.cend(),
                  dest_values.begin());
    }
    //=============================================================================
    void vec_axpy(real_t a, const Vector &x, Vector &y)
    {
        SFEM_CHECK_SIZES(x.block_size(), y.block_size());
        SFEM_CHECK_SIZES(x.n_local(), y.n_local());
        for (int i = 0; i < x.n_local(); i++)
        {
            for (int j = 0; j < x.block_size(); j++)
            {
                y(i, j) += a * x(i, j);
            }
        }
    }
}