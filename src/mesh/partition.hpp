#pragma once

#include "cell.hpp"
#include "../graph/partition.hpp"

namespace sfem::mesh
{
    /// @brief The criterion based on which two cells are
    /// considered to be linked/connected. For example,
    /// if shared_node is selected, two cells are considered
    /// linked if they share a single node.
    enum class PartitionCriterion
    {
        shared_facet,
        shared_node
    };

    /// @brief Create the cell partition
    /// @param cells The cells
    /// @param cell_to_node The cell-to-node connectivity
    /// @param partition_criteration The criterion based on which cells are
    /// considered connected
    /// @param partitioner_type The partitioner to be used
    /// @return The cell partition (as an indexmap)
    std::shared_ptr<IndexMap>
    create_cell_partition(const std::vector<Cell> &cells,
                          const graph::Connectivity &cell_to_node,
                          PartitionCriterion partition_criteration,
                          graph::partition::PartitionerType partitioner_type);

    /// @brief Partition a certain set of mesh entities using the existing cell partition.
    /// The resulting partition is constructed by the simple rule that out of all the cells
    /// connected to an entity, the cell that owns it the entity is the one with the highest
    /// global index. By extension, the entity is owned by that cell's owner process.
    /// @param cell_im The cell indexmap (cell partition)
    /// @param cell_to_entity The (old) cell-to-entity connectivity
    /// @return The entity partition and the (new) renumbered cell-to-entity
    /// connectivity.
    std::pair<std::shared_ptr<IndexMap>,
              std::shared_ptr<graph::Connectivity>>
    create_entity_partition(const IndexMap &cell_im,
                            const graph::Connectivity &cell_to_entity);
}