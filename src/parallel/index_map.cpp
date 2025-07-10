#include "index_map.hpp"
#include "mpi.hpp"
#include "../base/error.hpp"
#include <ranges>
#include <algorithm>
#include <numeric>

namespace sfem
{
    //=============================================================================
    IndexMap::IndexMap(int n_owned)
        : local_to_global_(n_owned),
          ghost_owners_(),
          global_to_local_()

    {
        for (int i = 0; i < n_owned; i++)
        {
            local_to_global_[i] = i;
            global_to_local_[i] = i;
        }
    }
    //=============================================================================
    IndexMap::IndexMap(std::vector<int> &&global_idxs,
                       std::vector<int> &&ghost_owners)
        : local_to_global_(std::move(global_idxs)),
          ghost_owners_(std::move(ghost_owners))
    {
        // Create the global-to-local mapping
        for (int i = 0; i < n_local(); i++)
        {
            global_to_local_[local_to_global_[i]] = i;
        }
    }
    //=============================================================================
    int IndexMap::n_owned() const
    {
        return static_cast<int>(local_to_global_.size() - ghost_owners_.size());
    }
    //=============================================================================
    int IndexMap::n_ghost() const
    {
        return static_cast<int>(ghost_owners_.size());
    }
    //=============================================================================
    int IndexMap::n_local() const
    {
        return static_cast<int>(local_to_global_.size());
    }
    //=============================================================================
    int IndexMap::n_global() const
    {
        if (mpi::n_procs() == 1)
        {
            return n_owned();
        }
        else
        {
            return mpi::reduce(n_owned(), mpi::ReduceOperation::sum);
        }
    }
    //=============================================================================
    std::vector<int> IndexMap::owned_idxs() const
    {
        return {local_to_global_.cbegin(),
                local_to_global_.cbegin() + n_owned()};
    }
    //=============================================================================
    std::vector<int> IndexMap::ghost_idxs() const
    {
        return {local_to_global_.cbegin() + n_owned(),
                local_to_global_.cend()};
    }
    //=============================================================================
    std::vector<int> IndexMap::local_idxs() const
    {
        return local_to_global_;
    }
    //=============================================================================
    std::vector<int> IndexMap::ghost_owners() const
    {
        return ghost_owners_;
    }
    //=============================================================================
    int IndexMap::local_to_global(int local_idx) const
    {
        SFEM_CHECK_INDEX(local_idx, n_local());
        return local_to_global_[local_idx];
    }
    //=============================================================================
    std::vector<int> IndexMap::local_to_global(const std::span<const int> local_idxs) const
    {
        std::vector<int> global_idxs(local_idxs.size());
        for (std::size_t i = 0; i < local_idxs.size(); i++)
        {
            global_idxs[i] = local_to_global(local_idxs[i]);
        }
        return global_idxs;
    }
    //=============================================================================
    int IndexMap::global_to_local(int global_idx) const
    {
        if (global_to_local_.contains(global_idx))
        {
            return global_to_local_.at(global_idx);
        }
        else
        {
            return -1;
        }
    }
    //=============================================================================
    std::vector<int> IndexMap::global_to_local(const std::span<const int> global_idxs) const
    {
        std::vector<int> local_idxs(global_idxs.size());
        for (std::size_t i = 0; i < global_idxs.size(); i++)
        {
            local_idxs[i] = global_to_local(global_idxs[i]);
        }
        return local_idxs;
    }
    //=============================================================================
    int IndexMap::get_owner(int local_idx) const
    {
        if (local_idx < n_owned())
        {
            return mpi::rank();
        }
        else if (local_idx < n_local())
        {
            return ghost_owners_[local_idx - n_owned()];
        }
        else
        {
            SFEM_CHECK_INDEX(local_idx, n_local());
            return -1;
        }
    }
    //=============================================================================
    bool IndexMap::is_ghost(int local_idx) const
    {
        if (get_owner(local_idx) != mpi::rank())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    //=============================================================================
    IndexMap IndexMap::renumber() const
    {
        // Get rank and number of process
        int proc_rank = mpi::rank();
        int n_procs = mpi::n_procs();

        // If serial, return early
        if (n_procs == 1)
        {
            return IndexMap(n_owned());
        }

        // Vector of renumbered indices present in this process
        // The first n_owned indices are owned by this process,
        // while the last n_ghost are owned by others
        std::vector<int> global_idxs(n_owned() + n_ghost());

        // Owned
        {
            // Compute the displacement for each process
            std::vector<int> send_buffer(n_procs, n_owned());
            std::vector<int> owners(n_procs);
            std::ranges::iota(owners, 0);
            auto [recv_buffer, _, __] = mpi::send_to_dest<int>(send_buffer, owners);

            // Compute the renumbered owned indices
            int disp = std::accumulate(recv_buffer.cbegin(), recv_buffer.cbegin() + proc_rank, 0);
            std::iota(global_idxs.begin(), global_idxs.begin() + n_owned(), disp);
        }

        // Ghost
        {
            // Send the ghost indices to their owner processes
            std::span<const int> ghosts(local_to_global_.cbegin() + n_owned(), n_ghost());
            auto [recv_buffer, recv_counts, recv_displs] = mpi::send_to_dest<int>(ghosts, ghost_owners_);

            // Renumber the received ghost indices
            for (std::size_t i = 0; i < recv_buffer.size(); i++)
            {
                int local_idx = global_to_local_.at(recv_buffer[i]);
                recv_buffer[i] = global_idxs[local_idx];
            }

            // Send the renumbered ghost indices back to the process that sent them
            std::vector<int> owners;
            for (int i = 0; i < n_procs; i++)
            {
                for (int j = 0; j < recv_counts[i]; j++)
                {
                    owners.emplace_back(i);
                }
            }
            std::tie(recv_buffer, recv_counts, recv_displs) = mpi::send_to_dest<int>(recv_buffer, owners);

            // Store the renumbered ghost indices
            std::fill(recv_counts.begin(), recv_counts.end(), 0);
            for (int i = n_owned(); i < n_owned() + n_ghost(); i++)
            {
                int owner = ghost_owners_[i - n_owned()];
                global_idxs[i] = recv_buffer[recv_displs[owner] + recv_counts[owner]++];
            }
        }

        return IndexMap(std::move(global_idxs), ghost_owners());
    }
}