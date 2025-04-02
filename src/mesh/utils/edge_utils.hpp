#pragma once

#include "../cell.hpp"
#include "../../graph/connectivity.hpp"
#include "../../parallel/index_map.hpp"
#include <map>
#include <memory>

namespace sfem::mesh::utils
{
    class EdgeMap
    {
    public:
        /// @brief Insert an edge into the map, if it's not already included
        /// @param edge_nodes The nodes of the edge. They need not be in order
        /// @return The edge's index
        int insert(std::array<int, 2> edge_nodes);

    private:
        std::map<std::array<int, 2>, int> map_;
    };

    /// @brief Extract the cell edges
    /// @param cells The cells
    /// @param cell_to_node The cell-to-node connectivity
    /// @return The cell-to-edge connectivity
    std::shared_ptr<graph::Connectivity>
    extract_edges(const std::vector<Cell> &cells,
                  const graph::Connectivity &cell_to_node);

    /// @brief Generate the edge-to-node connectivity
    /// @param cells The cells
    /// @param cell_to_edge The cell-to-edge connectivity
    /// @param cell_to_node The cell-to-node connectivity
    /// @param cell_index_map The cell index map
    /// @return The edge-to-node connectivity
    /// @note The edge's orientation (i.e. the ordering of its nodes)
    /// matches that of the adjacent cell with the highest global index
    std::shared_ptr<graph::Connectivity>
    edge_to_node(const std::vector<Cell> &cells,
                 const graph::Connectivity &cell_to_edge,
                 const graph::Connectivity &cell_to_node,
                 const IndexMap &cell_index_map);
}