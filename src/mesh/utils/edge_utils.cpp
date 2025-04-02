#include "edge_utils.hpp"
#include "../../base/timer.hpp"
#include "../../base/error.hpp"
#include <numeric>

namespace sfem::mesh::utils
{
    //=========================================================================
    int EdgeMap::insert(std::array<int, 2> edge_nodes)
    {
        // Sort the edge nodes in ascending order,
        // if necessary
        if (edge_nodes[0] > edge_nodes[1])
        {
            std::swap(edge_nodes[1], edge_nodes[0]);
        }

        // Insert only if the edge is not already registered
        // In either case, return the edge's index
        int edge_idx;
        if (!map_.contains(edge_nodes))
        {
            edge_idx = static_cast<int>(map_.size());
            map_.insert(std::make_pair(edge_nodes, edge_idx));
        }
        else
        {
            edge_idx = map_[edge_nodes];
        }

        return edge_idx;
    }
    //=========================================================================
    std::shared_ptr<graph::Connectivity>
    extract_edges(const std::vector<Cell> &cells,
                  const graph::Connectivity &cell_to_node)
    {
        Timer timer;

        SFEM_CHECK_SIZES(cells.size(), cell_to_node.n_primary());

        // Compute the cell-to-edge connectivity offsets
        std::vector<int> cell_n_edges(cells.size(), 0);
        for (std::size_t i = 0; i < cells.size(); i++)
        {
            cell_n_edges[i] = cell_num_edges(cells[i].type);
        }
        std::vector<int> cell_edge_offsets(cells.size() + 1, 0);
        std::inclusive_scan(cell_n_edges.cbegin(),
                            cell_n_edges.cend(),
                            cell_edge_offsets.begin() + 1);

        // Fill the cell-to-edge connectivity array
        std::vector<int> cell_edge_array(cell_edge_offsets.back(), 0);
        std::fill(cell_n_edges.begin(), cell_n_edges.end(), 0);
        EdgeMap edge_map;
        for (int i = 0; i < cell_to_node.n_primary(); i++)
        {
            auto cell = cells[i];
            auto cell_nodes = cell_to_node.links(i);
            for (int j = 0; j < cell_num_edges(cell.type); j++)
            {
                auto edge_ordering = cell_edge_ordering(cell.type, j);
                std::array<int, 2> edge_nodes;
                for (int k = 0; k < 2; k++)
                {
                    edge_nodes[k] = cell_nodes[edge_ordering[k]];
                }
                cell_edge_array[cell_edge_offsets[i] + cell_n_edges[i]++] = edge_map.insert(edge_nodes);
            }
        }

        return std::make_shared<graph::Connectivity>(std::move(cell_edge_offsets),
                                                     std::move(cell_edge_array));
    }
    //=========================================================================
    std::shared_ptr<graph::Connectivity>
    edge_to_node(const std::vector<Cell> &cells,
                 const graph::Connectivity &cell_to_edge,
                 const graph::Connectivity &cell_to_node,
                 const IndexMap &cell_index_map)
    {
        SFEM_CHECK_SIZES(cells.size(), cell_to_edge.n_primary());
        SFEM_CHECK_SIZES(cells.size(), cell_to_node.n_primary());
        SFEM_CHECK_SIZES(cells.size(), cell_index_map.n_local());

        // Compute the edge-to-node connectivity offsets
        std::vector<int> edge_node_offsets(cell_to_edge.n_secondary() + 1, 2);
        std::exclusive_scan(edge_node_offsets.cbegin(),
                            edge_node_offsets.cend(),
                            edge_node_offsets.begin(), 0);

        // Fill the edge-to-node connectivity array
        std::vector<int> edge_node_array(edge_node_offsets.back(), 0);
        auto edge_to_cell = cell_to_edge.invert();
        for (int i = 0; i < edge_to_cell.n_primary(); i++)
        {
            auto edge_cells = cell_index_map.local_to_global(edge_to_cell.links(i));
            auto owner_cell_idx = cell_index_map.global_to_local(std::ranges::max(edge_cells));
            auto owner_cell_type = cells[owner_cell_idx].type;
            auto owner_cell_nodes = cell_to_node.links(owner_cell_idx);
            auto edge_rel_idx = cell_to_edge.relative_index(owner_cell_idx, i);
            auto edge_ordering = cell_edge_ordering(owner_cell_type, edge_rel_idx);
            for (int j = 0; j < 2; j++)
            {
                edge_node_array[edge_node_offsets[i] + j] = owner_cell_nodes[edge_ordering[j]];
            }
        }

        return std::make_shared<graph::Connectivity>(std::move(edge_node_offsets),
                                                     std::move(edge_node_array));
    }
}