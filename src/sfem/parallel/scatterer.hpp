#pragma once

#include <sfem/parallel/index_map.hpp>
#include <sfem/parallel/mpi.hpp>
#include <sfem/base/error.hpp>
#include <functional>
#include <memory>

namespace sfem
{
    /// @brief This class facilitates communication of ghost index values
    /// for a given index map.
    ///
    /// A forward scatter consists of each process sending values for locally owned
    /// indices which are ghosted on other processes. These indices are computed and
    /// stored during construction, to avoid repeated communication overhead. It is
    /// important to note than more than one processes might ghost the same locally
    /// owned index. Due to this, a forward scatter index might have more than one
    /// destination process.
    ///
    /// A reverse scatter consists of each process sending its stored ghost index
    /// values to their respective owner processes.
    template <typename T>
    class Scatterer
    {
    public:
        /// @brief Create a scatterer for a given index map
        /// @param index_map
        Scatterer(std::shared_ptr<const IndexMap> index_map) : index_map_(index_map)
        {
            // Compute the forward scatter indices
            auto [fwd_idxs_global, counts, displs] = (mpi::send_to_dest<int>(index_map_->ghost_idxs(),
                                                                             index_map_->ghost_owners()));
            fwd_idxs_ = index_map_->global_to_local(fwd_idxs_global);

            // Store the destination process for each forward scatter index
            fwd_dest_.resize(fwd_idxs_.size());
            for (int i = 0; i < mpi::n_procs(); i++)
            {
                for (int j = 0; j < counts[i]; j++)
                {
                    fwd_dest_[displs[i] + j] = i;
                }
            }

            // Compute the reverse scatter indices
            auto rev_idxs_global = std::get<0>(mpi::send_to_dest<int>(fwd_idxs_global, fwd_dest_));
            rev_idxs_ = index_map_->global_to_local(rev_idxs_global);
        }

        std::shared_ptr<const IndexMap> index_map() const
        {
            return index_map_;
        }

        std::vector<int> forward_idxs() const
        {
            return fwd_idxs_;
        }

        std::vector<int> reverse_idxs() const
        {
            return rev_idxs_;
        }

        /// @brief Send values for locally owned indices, while receiving values for ghost indices
        /// @param values Buffer containing values for both locally owned and ghost indices
        /// @param bs Block size
        /// @param op Binary operation to modify values
        void forward(std::span<T> values, int bs, std::function<void(T &, T)> op) const
        {
            SFEM_CHECK_SIZES(values.size(), index_map_->n_local() * bs);

            // Create the send buffer, i.e. store the values for the
            // forward scatter indices in the correct order
            std::vector<T> send_buffer(fwd_idxs_.size() * bs);
            for (std::size_t i = 0; i < fwd_idxs_.size(); i++)
            {
                for (int j = 0; j < bs; j++)
                {
                    send_buffer[i * bs + j] = values[fwd_idxs_[i] * bs + j];
                }
            }

            // Send the values for the forward scatter indices.
            // Receive values for ghost indices
            auto recv_buffer = std::get<0>(mpi::send_to_dest<T>(send_buffer, fwd_dest_, bs));

            // Modify ghost index values
            for (std::size_t i = 0; i < rev_idxs_.size(); i++)
            {
                for (int j = 0; j < bs; j++)
                {
                    op(values[rev_idxs_[i] * bs + j], recv_buffer[i * bs + j]);
                }
            }
        }

        /// @brief Send values for ghosted indices, while receiving values for locally owned indices
        /// @param values Buffer containing values for both locally owned and ghost indices
        /// @param bs Block size
        /// @param op Binary operation to modify values
        void reverse(std::span<T> values, int bs, std::function<void(T &, T)> op) const
        {
            SFEM_CHECK_SIZES(values.size(), index_map_->n_local() * bs);

            // Quick access
            const int n_owned = index_map_->n_owned();

            // Create the send buffer (ghost index values)
            std::span<T> send_buffer = {values.begin() + n_owned * bs, values.end()};

            auto recv_buffer = std::get<0>(mpi::send_to_dest<T>(send_buffer,
                                                                index_map_->ghost_owners(), bs));

            // Modify locally owned index values
            for (std::size_t i = 0; i < fwd_idxs_.size(); i++)
            {
                for (int j = 0; j < bs; j++)
                {

                    op(values[fwd_idxs_[i] * bs + j], recv_buffer[i * bs + j]);
                }
            }
        }

    private:
        /// @brief Index map
        std::shared_ptr<const IndexMap> index_map_;

        /// @brief Forward scatter indices, i.e. locally owned
        /// indices that are ghosted on other processes
        std::vector<int> fwd_idxs_;

        /// @brief Destination process for each forward scatter index
        std::vector<int> fwd_dest_;

        /// @brief Reverse scatter indices, i.e. ghosted indices
        /// @note They are stored in a different order compared
        /// to the index map.
        std::vector<int> rev_idxs_;
    };
}