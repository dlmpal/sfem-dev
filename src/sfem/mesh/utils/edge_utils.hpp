#pragma once

#include <sfem/mesh/cell.hpp>
#include <sfem/graph/connectivity.hpp>
#include <sfem/parallel/index_map.hpp>
#include <map>
#include <memory>

namespace sfem::mesh::utils
{
    class EdgeMap
    {
    public:
        /// @brief Insert an edge into the map, if it's not already included
        /// @param edge_nodes Edge nodes. They need not be in order
        /// @return The edge's index
        int insert(std::array<int, 2> edge_nodes);

        /// @brief Insert an edge into the map, with a specified index
        /// @param edge_nodes Edge nodes. They need not be in order
        /// @param edge_idx The edge's index
        void insert(std::array<int, 2> edge_nodes, int edge_idx);

        /// @brief Get the the index of an edge, given its nodes
        /// @param edge_nodes Edge nodes. They need not be in order
        /// @return The edge index
        int at(std::array<int, 2> edge_nodes);

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

    /// @brief Generate the face-to-edge connectivity
    /// from existing connectivities
    /// @param cells The cells
    /// @param faces The faces
    /// @param cell_to_edge The cell-to-edge connectivity
    /// @param cell_to_node The cell-to-node connectivity
    /// @param face_to_node The face-to-node connectivity
    /// @return The face-to-edge connectivity
    std::shared_ptr<graph::Connectivity>
    face_to_edge(const std::vector<Cell> &cells,
                 const std::vector<Cell> &faces,
                 const graph::Connectivity &cell_to_edge,
                 const graph::Connectivity &cell_to_node,
                 const graph::Connectivity &face_to_node);
}