#include "partition.hpp"
#include <sfem/mesh/topology.hpp>
#include <sfem/mesh/utils/face_utils.hpp>
#include <sfem/mesh/utils/edge_utils.hpp>
#include <sfem/parallel/mpi.hpp>
#include <sfem/base/error.hpp>
#include <numeric>
#include <ranges>

namespace sfem::mesh
{
    //=============================================================================
    std::shared_ptr<IndexMap> create_cell_partition(const std::vector<Cell> &cells,
                                                    const graph::Connectivity &cell_to_node,
                                                    PartitionCriterion partition_criterion,
                                                    graph::partition::PartitionerType partitioner_type)
    {
        SFEM_CHECK_SIZES(cells.size(), cell_to_node.n_primary());

        // The cell partition
        std::vector<int> cell_partition;

        // In order for the partition to be computed,
        // the cell-to-cell connectivity must first be extracted.
        // The choice of partition criterion results in different cell-to-cell
        // connectivities, and this partitions. This work,
        // as well as the partitioning, are performed only by the root process
        graph::Connectivity cell_to_cell;
        if (mpi::rank() == mpi::root())
        {
            int dim = cell_dim(cells.front().type);

            if (partition_criterion == PartitionCriterion::shared_facet and dim == 3)
            {
                auto [_, cell_to_face] = utils::extract_faces(cells, cell_to_node);
                cell_to_cell = cell_to_face->primary_to_primary();
            }
            if (partition_criterion == PartitionCriterion::shared_facet and dim == 2)
            {
                auto cell_to_edge = utils::extract_edges(cells, cell_to_node);
                cell_to_cell = cell_to_edge->primary_to_primary();
            }
            else
            {
                cell_to_cell = cell_to_node.primary_to_primary();
            }

            cell_partition = graph::partition::create_partition(cell_to_cell,
                                                                mpi::n_procs(),
                                                                partitioner_type);
        }

        return graph::partition::distribute_partition(cell_to_cell, cell_partition);
    }
    //=========================================================================
    std::pair<std::shared_ptr<IndexMap>,
              std::shared_ptr<graph::Connectivity>>
    create_entity_partition(const IndexMap &cell_im,
                            const graph::Connectivity &cell_to_entity)
    {
        // Check that sizes match
        SFEM_CHECK_SIZES(cell_im.n_local(), cell_to_entity.n_primary());

        // Get the inverse connectivity
        auto entity_to_cell = cell_to_entity.invert();

        // Renumber the entities locally,
        // so that they are ordered as follows:
        // [0, owned_1, owned_2, ..., n_owned_local-1,
        // ghost_1, ghost_2, ..., n_owned_local+n_ghost)
        // Thus, this vector maps an entity's old local index
        // to the renumbered one
        std::vector<int> local_renumbering(entity_to_cell.n_primary());
        std::vector<bool> is_ghost(entity_to_cell.n_primary(), false);
        int n_owned = 0;
        int n_ghost = 0;

        // This struct contains the three following vectors:
        // 1) The global owner cell index for each ghost entity
        // 2) The ghost entity's index in its owner cell
        // 3) The owner process for the owner cell
        struct GhostEntities
        {
            std::vector<int> owner_cell_idxs;
            std::vector<int> idx_in_owner_cell;
            std::vector<int> owner_procs;
        };

        // Renumber all entities locally
        // For the ghost entities, store the relevant info
        GhostEntities ghosts;
        for (int i = 0; i < entity_to_cell.n_primary(); i++)
        {
            // Get all cells linked to the entity
            auto linked_cells_global_idxs = cell_im.local_to_global(entity_to_cell.links(i));

            // Get the owner cell for the entity (both in local and global indexing)
            // Out of all cells linked to the entity, the owner is that with the
            // highest global index
            int owner_cell_global_idx = std::ranges::max(linked_cells_global_idxs);
            int owner_cell_local_idx = cell_im.global_to_local(owner_cell_global_idx);

            // Get the relative index of the entity in the owner cell
            int idx_in_owner_cell = cell_to_entity.relative_index(owner_cell_local_idx, i);

            // Get the owner process for the entity's owner cell
            int owner_proc = cell_im.get_owner(owner_cell_local_idx);

            if (owner_proc == mpi::rank())
            {
                local_renumbering[i] = n_owned++;
            }
            else
            {
                local_renumbering[i] = n_ghost++;
                is_ghost[i] = true;
                ghosts.owner_cell_idxs.emplace_back(owner_cell_global_idx);
                ghosts.idx_in_owner_cell.emplace_back(idx_in_owner_cell);
                ghosts.owner_procs.emplace_back(owner_proc);
            }
        }

        // Correctly offset the ghost entities
        for (int i = 0; i < entity_to_cell.n_primary(); i++)
        {
            if (is_ghost[i])
            {
                local_renumbering[i] += n_owned;
            }
        }

        // Apply the renumbering to the cell-to-entity connectivity
        auto cell_entity_array = cell_to_entity.array();
        for (std::size_t i = 0; i < cell_entity_array.size(); i++)
        {
            cell_entity_array[i] = local_renumbering[cell_entity_array[i]];
        }
        auto cell_to_entity_re = std::make_shared<graph::Connectivity>(cell_to_entity.offsets(),
                                                                       std::move(cell_entity_array));
        auto entity_to_cell_re = cell_to_entity_re->invert();

        // Now that the entities are correctly numbered locally,
        // the process of building the indexmap can begin.
        // First, the global offset for this process has to be
        // computed, so that the owned entities of this process
        // are in the following global range:
        // [offset, n_owned + offset)
        std::vector<int> send_buf(mpi::n_procs(), n_owned);
        std::vector<int> send_dest(mpi::n_procs(), 0);
        std::iota(send_dest.begin(), send_dest.end(), 0);
        auto recv_buf = std::get<0>(mpi::send_to_dest<int>(send_buf, send_dest));

        int offset = std::accumulate(recv_buf.cbegin(),
                                     recv_buf.cbegin() + mpi::rank(), 0);

        // Start building the local-to-global mapping for the entity indexmap
        std::vector<int> local_to_global(n_owned + n_ghost, -1);
        std::iota(local_to_global.begin(), local_to_global.begin() + n_owned, offset);

        // The correct numbering of the ghost entities must be communicated to this process
        // by their owner processes. In order for this to be achieved, each process sends
        // (for each ghost entity) the owner cell's global index and the index of the entity
        // in that cell, to the owner cell's process.
        auto [recv_ghost_entity_owner_cell_idxs,
              recv_ghost_entity_counts,
              recv_ghost_entity_displs] = mpi::send_to_dest<int>(ghosts.owner_cell_idxs,
                                                                 ghosts.owner_procs);

        auto [recv_ghost_entity_idx_in_owner_cell,
              _0, _1] = mpi::send_to_dest<int>(ghosts.idx_in_owner_cell,
                                               ghosts.owner_procs);

        // Indirect ghosts are ghost entities to which a process has assigned
        // an incorrect owner, due to the nature of the cell partitioning.
        // In order for this to be resolved, the current process,
        // which has received the indirect ghosts, must ask the actual
        // owner process of these entities for their correct global numbering.
        // After this is communicated to the current process, the correct global
        // indices of all ghost entities it received (indirect included) are sent
        // back the process which sent them, along with an array indicating the *actual*
        // owner process for each ghost.
        GhostEntities indirect_ghosts;
        std::vector<int> indirect_ghost_pos;
        std::vector<int> send_ghost_entity_idxs(recv_ghost_entity_owner_cell_idxs.size());
        std::vector<int> send_ghost_entity_owner_proc(recv_ghost_entity_owner_cell_idxs.size());
        std::vector<int> send_ghost_entity_dest(recv_ghost_entity_owner_cell_idxs.size());
        for (int i = 0; i < mpi::n_procs(); i++)
        {
            for (int j = 0; j < recv_ghost_entity_counts[i]; j++)
            {
                // The position of the entity in the received buffer
                int pos = recv_ghost_entity_displs[i] + j;

                // The local index of the entity's perceived owner cell
                int owner_cell_local_idx = cell_im.global_to_local(recv_ghost_entity_owner_cell_idxs[pos]);

                // The entity's relative index in the perceived owner cell
                int idx_in_owner_cell = recv_ghost_entity_idx_in_owner_cell[pos];

                // The entity's local index
                int entity_local_idx = cell_to_entity_re->links(owner_cell_local_idx)[idx_in_owner_cell];

                // Get the local and global index of the actual owner cell
                auto linked_cells = cell_im.local_to_global(entity_to_cell_re.links(entity_local_idx));
                int owner_cell_global_idx = std::ranges::max(linked_cells);
                owner_cell_local_idx = cell_im.global_to_local(owner_cell_global_idx);

                // Get the relative index of the entity in the actual owner cell
                idx_in_owner_cell = cell_to_entity_re->relative_index(owner_cell_local_idx, entity_local_idx);

                // Store the entity's global index and its dest. process
                if (cell_im.get_owner(owner_cell_local_idx) == mpi::rank())
                {
                    send_ghost_entity_idxs[pos] = local_to_global[entity_local_idx];
                }
                else
                {
                    indirect_ghosts.owner_cell_idxs.emplace_back(owner_cell_global_idx);
                    indirect_ghosts.idx_in_owner_cell.emplace_back(idx_in_owner_cell);
                    indirect_ghosts.owner_procs.emplace_back(cell_im.get_owner(owner_cell_local_idx));
                    indirect_ghost_pos.emplace_back(pos);
                }
                send_ghost_entity_owner_proc[pos] = cell_im.get_owner(owner_cell_local_idx);
                send_ghost_entity_dest[pos] = i;
            }
        }

        auto [recv_indirect_ghost_entity_owner_cell_idxs,
              recv_indirect_ghost_entity_counts,
              recv_indirect_ghost_entity_displs] = mpi::send_to_dest<int>(indirect_ghosts.owner_cell_idxs,
                                                                          indirect_ghosts.owner_procs);
        auto [recv_indirect_ghost_entity_idx_in_owner_cell,
              _2, _3] = mpi::send_to_dest<int>(indirect_ghosts.idx_in_owner_cell,
                                               indirect_ghosts.owner_procs);

        // Each process assigns the correct global index to each indirect ghost it has received
        std::vector<int> send_indirect_ghost_entity_idxs(recv_indirect_ghost_entity_owner_cell_idxs.size());
        std::vector<int> send_indirect_ghost_entity_dest(recv_indirect_ghost_entity_owner_cell_idxs.size());
        for (int i = 0; i < mpi::n_procs(); i++)
        {
            for (int j = 0; j < recv_indirect_ghost_entity_counts[i]; j++)
            {
                // The position of the entity in the received buffer
                int pos = recv_indirect_ghost_entity_displs[i] + j;

                // The local index of the entity's owner cell
                int owner_cell_global_idx = recv_indirect_ghost_entity_owner_cell_idxs[pos];
                int owner_cell_local_idx = cell_im.global_to_local(owner_cell_global_idx);

                // The entity's index in the owner cell
                int idx_in_owner_cell = recv_indirect_ghost_entity_idx_in_owner_cell[pos];

                // The entity's local index
                int entity_local_idx = cell_to_entity_re->links(owner_cell_local_idx)[idx_in_owner_cell];

                // Store the entity's global index and its dest. process
                send_indirect_ghost_entity_idxs[pos] = local_to_global[entity_local_idx];
                send_indirect_ghost_entity_dest[pos] = i;
            }
        }

        // After the indirect ghost entities are assigned correct global indices from their actual owner process,
        // their are inserted into the send buffer, along with all other ghost entities that have been
        // received by the current process
        std::vector<int> recv_indirect_ghost_entity_idxs;
        std::tie(recv_indirect_ghost_entity_idxs,
                 recv_indirect_ghost_entity_counts,
                 recv_indirect_ghost_entity_displs) = mpi::send_to_dest<int>(send_indirect_ghost_entity_idxs,
                                                                             send_indirect_ghost_entity_dest);
        for (std::size_t i = 0; i < indirect_ghost_pos.size(); i++)
        {
            int owner_proc = indirect_ghosts.owner_procs[i];
            int pos = recv_indirect_ghost_entity_displs[owner_proc]++;
            int global_idx = recv_indirect_ghost_entity_idxs[pos];
            send_ghost_entity_idxs[indirect_ghost_pos[i]] = global_idx;
        }

        // The correct global indices for the ghost entities and their actual owner processes
        // are sent back to their destination processes.
        std::vector<int> ghost_entity_global_idxs;
        std::vector<int> ghost_entity_actual_owner_procs;
        std::tie(ghost_entity_global_idxs,
                 recv_ghost_entity_counts,
                 recv_ghost_entity_displs) = mpi::send_to_dest<int>(send_ghost_entity_idxs, send_ghost_entity_dest);
        std::tie(ghost_entity_actual_owner_procs,
                 _0, _1) = mpi::send_to_dest<int>(send_ghost_entity_owner_proc, send_ghost_entity_dest);

        // Finally, finish the local-to-global mapping by
        // assigning the correct global index to each ghost entity.
        // Also, store the actual owner process for each ghost entity
        std::vector<int> ghost_owners(n_ghost);
        for (int i = 0; i < n_ghost; i++)
        {
            // The local index of the entity's owner cell
            int owner_cell_local_idx = cell_im.global_to_local(ghosts.owner_cell_idxs[i]);

            // The entity's index in the owner cell
            int idx_in_owner_cell = ghosts.idx_in_owner_cell[i];

            // The entity's (perceived) owner process
            int owner_proc = ghosts.owner_procs[i];

            // The position of the entity in the received buffer
            int pos = recv_ghost_entity_displs[owner_proc]++;

            // The local and global index of the entity
            int entity_local_idx = cell_to_entity_re->links(owner_cell_local_idx)[idx_in_owner_cell];
            int entity_global_idx = ghost_entity_global_idxs[pos];

            // Store the entity's global index and its owner process
            local_to_global[entity_local_idx] = entity_global_idx;
            ghost_owners[entity_local_idx - n_owned] = ghost_entity_actual_owner_procs[pos];
        }

        return {std::make_shared<IndexMap>(std::move(local_to_global),
                                           std::move(ghost_owners)),
                cell_to_entity_re};
    }
}