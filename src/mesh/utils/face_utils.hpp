#pragma once

#include "../cell.hpp"
#include "../../graph/connectivity.hpp"
#include "../../parallel/index_map.hpp"
#include <algorithm>
#include <map>
#include <memory>

namespace sfem::mesh::utils
{
    class FaceMap
    {
    public:
        /// @brief Insert a face into the map, if it's not already included.
        /// @param face_nodes The nodes of the face. They need not be in order
        /// @param face_type The face's type
        /// @return The face's index
        int insert(const std::span<const int> face_nodes, CellType face_type);

        /// @brief Get the face's type and its index, given its nodes.
        /// @param face_nodes The nodes of the face. They need not be in order
        /// @return The face type and its index
        /// @note Produces an error if the face is not found
        /// @note The size of face_nodes must equal the no. nodes for the face's cell type
        std::optional<std::pair<CellType, int>> at(const std::span<const int> face_nodes) const;

        /// @brief Get the faces
        std::vector<CellType> get_faces() const;

        // private:
        //  A face is uniquely identified by its nodes, sorted in
        //  ascending order. A face has either 3 (triangle) or 4 (quadrilateral) nodes.
        //  Thus an array of the max possible size (4), along with the
        //  actual number of nodes are stored. The nodes should be already
        //  sorted.
        struct _Key
        {
            std::array<int, 4> data{};
            int size{};

            friend bool operator<(const _Key &lhs, const _Key &rhs)
            {
                for (int i = 0; i < std::min(lhs.size, rhs.size); i++)
                {
                    if (lhs.data[i] < rhs.data[i])
                    {
                        return true;
                    }

                    if (lhs.data[i] > rhs.data[i])
                    {
                        return false;
                    }
                }

                return false;
            };
        };

        // For each face, its type and index are stored
        using _Value = std::pair<CellType, int>;

        /// @brief The actual face map
        std::map<_Key, _Value> map_;
    };

    /// @brief Extract the cell faces
    /// @param cells The cells
    /// @param cell_to_node The cell-to-node connectivity
    /// @return The face types and the cell-to-face connectivity
    std::pair<std::vector<CellType>,
              std::shared_ptr<graph::Connectivity>>
    extract_faces(const std::vector<Cell> &cells,
                  const graph::Connectivity &cell_to_node);

    /// @brief Generate the face-to-node connectivity
    /// @param cells The cells
    /// @param cell_to_face The cell-to-face connectivity
    /// @param cell_to_node The cell-to-node connectivity
    /// @param cell_index_map The cell index map
    /// @return The face-to-node connectivity
    /// @note The face's orientation (i.e. the ordering of its nodes)
    /// matches that of the adjacent cell with the highest global index
    std::shared_ptr<graph::Connectivity>
    face_to_node(const std::vector<Cell> &cells,
                 const graph::Connectivity &cell_to_face,
                 const graph::Connectivity &cell_to_node,
                 const IndexMap &cell_index_map);
}