#pragma once

#include <vector>
#include <unordered_map>
#include <span>

namespace sfem
{
    /// @brief This class is tasked with mapping a set of indices from indexing
    /// local to this process (or local indexing), to indexing used by all processes
    /// (or global indexing). Local indexing is the contiguous range [0, n_local),
    /// where n_local is the total number of indices present in this process, namely
    /// the sum of the number of indices owned by this process (or n_owned) and the
    /// number of indices present on this process but owned by others (or n_ghost).
    class IndexMap
    {
    public:
        /// @brief Create an IndexMap for serial execution
        /// @param n_owned Number of owned indices
        IndexMap(int n_owned = 0);

        /// @brief Create an IndexMap
        /// @param owned_idxs Range of global indices
        /// @param ghost_owners Owner process for each ghost index
        /// @note The first (global_idxs.size() - ghost_owners.size()) indices
        /// are considered owned, while the rest are considered ghosts
        IndexMap(std::vector<int> &&global_idxs,
                 std::vector<int> &&ghost_owners);

        /// @brief Get the number of owned indices
        int n_owned() const;

        /// @brief Get the number of ghost indices
        int n_ghost() const;

        /// @brief Get the number of local indices
        /// @note Sum of owned and ghost indices
        int n_local() const;

        /// @brief Get the global number of indices
        int n_global() const;

        /// @brief Get the owned indices
        std::vector<int> owned_idxs() const;

        /// @brief Get the ghost indices
        std::vector<int> ghost_idxs() const;

        /// @brief Get the local indices (owned + ghost)
        std::vector<int> local_idxs() const;

        /// @brief Get the owner process for each ghost index
        std::vector<int> ghost_owners() const;

        /// @brief Map an index from local to global indexing
        /// @note Raises an error if idx > n_local
        int local_to_global(int local_idx) const;

        /// @brief Map a range of indices from local to global indexing
        std::vector<int> local_to_global(const std::span<const int> local_idxs) const;

        /// @brief Map an index from global to local indexing
        /// @note Returns -1 if the index is not local
        int global_to_local(int global_idx) const;

        /// @brief Map a range of indices from global to local indexing
        std::vector<int> global_to_local(const std::span<const int> global_idxs) const;

        /// @brief Get the owner process for a (local) index
        int get_owner(int local_idx) const;

        /// @brief Check whether an index is a ghost for this process
        bool is_ghost(int local_idx) const;

        /// @brief Renumber the index map so that the owned indices of
        /// process 0 are in the range [0, n_owned_0), or more generally
        /// that the owned indices of process i are in the range
        /// [offset_i, offset_i + n_owned_i), where offset_i is the sum
        /// of n_owned_0 to n_owned_(i-1)
        /// @return The renumbered IndexMap
        IndexMap renumber() const;

    private:
        /// @brief Map from local to global indexing
        std::vector<int> local_to_global_;

        /// @brief Owner process for each ghost index
        std::vector<int> ghost_owners_;

        /// @brief Map from global to local indexing
        std::unordered_map<int, int> global_to_local_;
    };
}