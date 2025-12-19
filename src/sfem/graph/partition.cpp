#include "partition.hpp"
#include <sfem/parallel/mpi.hpp>
#include <sfem/base/error.hpp>
#include <sfem/base/timer.hpp>
#include <numeric>

#ifdef SFEM_HAS_METIS
#include <metis.h>
#endif

namespace sfem::graph::partition
{
    //=============================================================================
#ifdef SFEM_HAS_METIS
    std::vector<int> create_partition_metis(const graph::Connectivity &conn, int n_parts)
    {
        Timer timer;

        int n_vertices = conn.n_primary();
        int ncon = 1;
        int objval;
        std::vector<int> part(n_vertices, 0);
        int code;

        if (n_parts <= 8)
        {
            code = METIS_PartGraphRecursive(&n_vertices, &ncon,
                                            conn.offsets().data(), conn.array().data(),
                                            nullptr, nullptr, nullptr, &n_parts,
                                            nullptr, nullptr, nullptr,
                                            &objval, part.data());
        }
        else
        {
            code = METIS_PartGraphKway(&n_vertices, &ncon,
                                       conn.offsets().data(), conn.array().data(),
                                       nullptr, nullptr, nullptr, &n_parts,
                                       nullptr, nullptr, nullptr,
                                       &objval, part.data());
        }

        if (code != METIS_OK)
        {
            SFEM_ERROR(std::format("METIS exited with error {}\n", code));
        }

        return part;
    }
#endif
    //=============================================================================
    std::vector<int> create_partition(const Connectivity &conn, int n_parts,
                                      PartitionerType type)
    {
        // Return early for serial runs
        if (n_parts == 1)
        {
            return std::vector<int>(conn.n_primary(), 0);
        }

        auto type_to_string = [](PartitionerType type)
        {
            switch (type)
            {
            case PartitionerType::metis:
                return "METIS";
            default:
                SFEM_ERROR(std::format("Invalid partitioner type: {}\n", static_cast<int>(type)));
                return "";
            }
        };

        std::vector<int> owners;

#ifdef SFEM_HAS_METIS
        if (type == PartitionerType::metis)
        {
            owners = create_partition_metis(conn, n_parts);
        }
#endif

        // Perhaps not required
        // Keep just in case...
        if (owners.size() == 0)
        {
            SFEM_ERROR(std::format("Error creating partition with {} partitioner\n", type_to_string(type)));
        }

        return owners;
    }
    //=============================================================================
    std::shared_ptr<IndexMap> distribute_partition(const Connectivity &conn,
                                                   const std::span<const int> part)
    {
        // Owned indices
        std::vector<int> owned_idxs;

        // Owned indices destination process
        std::vector<int> owned_dest;

        // Ghost indices
        std::vector<int> ghost_idxs;

        // Ghost indices owners
        std::vector<int> ghost_owners;

        // Ghost indices destination process
        std::vector<int> ghost_dest;

        // The root process is tasked with filling the above vectors
        if (mpi::rank() == mpi::root())
        {
            SFEM_CHECK_SIZES(conn.n_primary(), part.size());

            // Get the number of partitions
            int n_parts = *std::max_element(part.cbegin(), part.cend()) + 1;

            // Vectors to store the number of owned and ghost indices
            // per partition
            std::vector<int> owned_count(n_parts, 0);
            std::vector<int> ghost_count(n_parts, 0);

            // Check if an index has already been included as a ghost
            // for a certain partition
            std::vector<std::unordered_map<int, bool>> is_included(n_parts);

            // Compute the number of owned and ghost indices for each partition
            for (int i = 0; i < conn.n_primary(); i++)
            {
                owned_count[part[i]]++;
                for (int j : conn.links(i))
                {
                    if (part[i] != part[j] and is_included[part[i]].contains(j) == false)
                    {
                        ghost_count[part[i]]++;
                        is_included[part[i]][j] = true;
                    }
                }
            }

            // Reset
            for (std::size_t i = 0; i < is_included.size(); i++)
            {
                is_included[i].clear();
            }

            // Create the displacement/offset vectors
            std::vector<int> owned_displ(n_parts, 0);
            std::vector<int> ghost_displ(n_parts, 0);
            std::exclusive_scan(owned_count.cbegin(), owned_count.cend(), owned_displ.begin(), 0);
            std::exclusive_scan(ghost_count.cbegin(), ghost_count.cend(), ghost_displ.begin(), 0);

            // Initialize the vectors to store the owned indices, the ghost indices and their owners
            // as well as the destination process for each index
            int n_owned_total = std::accumulate(owned_count.cbegin(), owned_count.cend(), 0);
            int n_ghost_total = std::accumulate(ghost_count.cbegin(), ghost_count.cend(), 0);
            owned_idxs.resize(n_owned_total);
            owned_dest.resize(n_owned_total);
            ghost_idxs.resize(n_ghost_total);
            ghost_owners.resize(n_ghost_total);
            ghost_dest.resize(n_ghost_total);

            // Reset the count vectors to use for indexing
            std::copy(owned_displ.cbegin(), owned_displ.cend(), owned_count.begin());
            std::copy(ghost_displ.cbegin(), ghost_displ.cend(), ghost_count.begin());

            // Fill the above vectors
            for (int i = 0; i < conn.n_primary(); i++)
            {
                owned_idxs[owned_count[part[i]]] = i;
                owned_dest[owned_count[part[i]]] = part[i];
                owned_count[part[i]]++;
                for (int j : conn.links(i))
                {
                    if (part[i] != part[j] and is_included[part[i]].contains(j) == false)
                    {
                        ghost_idxs[ghost_count[part[i]]] = j;
                        ghost_owners[ghost_count[part[i]]] = part[j];
                        ghost_dest[ghost_count[part[i]]] = part[i];
                        ghost_count[part[i]]++;
                        is_included[part[i]][j] = true;
                    }
                }
            }
        }

        // Broadcast data to all processes
        owned_idxs = mpi::distribute<int>(owned_idxs, owned_dest);
        ghost_idxs = mpi::distribute<int>(ghost_idxs, ghost_dest);
        ghost_owners = mpi::distribute<int>(ghost_owners, ghost_dest);

        // Append the ghost indices to the owned indices vector
        owned_idxs.insert(owned_idxs.end(), ghost_idxs.begin(), ghost_idxs.end());

        return std::make_shared<IndexMap>(std::move(owned_idxs), std::move(ghost_owners));
    }
}