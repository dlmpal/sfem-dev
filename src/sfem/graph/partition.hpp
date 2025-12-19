#pragma once

#include <sfem/graph/connectivity.hpp>
#include <sfem/parallel/index_map.hpp>
#include <memory>

/// @brief Functionality related to graph partitioning
namespace sfem::graph::partition
{
    enum class PartitionerType
    {
        metis
    };

    /// @brief Partition the given graph connectivity
    /// @note Should be called only by the root process
    std::vector<int> create_partition(const Connectivity &conn, int n_parts,
                                      PartitionerType type = PartitionerType::metis);

    /// @brief Distribute the connectivity partition data to all processes
    std::shared_ptr<IndexMap> distribute_partition(const Connectivity &conn,
                                                   const std::span<const int> part);
}