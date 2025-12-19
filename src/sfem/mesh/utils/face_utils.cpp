#include "face_utils.hpp"
#include <sfem/base/timer.hpp>
#include <sfem/base/error.hpp>
#include <numeric>

namespace sfem::mesh::utils
{
    //=========================================================================
    static void sort_4(std::span<int> data)
    {
        // Sort a range of either 2, 3, or 4 ints
        if (data.size() == 2)
        {
            if (data[0] > data[1])
            {
                std::swap(data[0], data[1]);
            }
        }
        else if (data.size() == 3)
        {
            if (data[0] > data[1])
            {
                std::swap(data[0], data[1]);
            }
            if (data[1] > data[2])
            {
                std::swap(data[1], data[2]);
            }
            if (data[0] > data[1])
            {
                std::swap(data[0], data[1]);
            }
        }
        else if (data.size() == 4)
        {
            sort_4({data.begin(), data.begin() + 3});
            if (data[3] < data[0])
            {
                std::rotate(data.rbegin(), data.rbegin() + 1, data.rend());
            }
            else if (data[3] < data[1])
            {
                std::swap(data[1], data[3]);
                std::swap(data[2], data[3]);
            }
            else if (data[3] < data[2])
            {
                std::swap(data[2], data[3]);
            }
        }
    }
    //=========================================================================
    int FaceMap::insert(const std::span<const int> face_nodes, CellType face_type)
    {
        SFEM_CHECK_SIZES(face_nodes.size(), cell_num_nodes(face_type));

        // Create the key
        _Key key;

        // Get no. nodes for the face
        key.size = cell_num_nodes(face_type);

        // Sort the nodes
        std::copy(face_nodes.cbegin(),
                  face_nodes.cbegin() + key.size,
                  key.data.begin());
        sort_4({key.data.begin(),
                key.data.begin() + key.size});

        // Insert only if the face is not already registered
        // Else change the face's region tag to 0,
        // indicating that the face is internal
        // In either case, return the face's index
        int face_idx;
        if (map_.contains(key) == false)
        {
            face_idx = static_cast<int>(map_.size());
            map_.insert(std::make_pair(key, std::make_pair(face_type, face_idx)));
        }
        else
        {
            face_idx = map_[key].second;
        }

        return face_idx;
    }
    //=========================================================================
    std::optional<std::pair<CellType, int>>
    FaceMap::at(const std::span<const int> face_nodes) const
    {
        // Create the key
        _Key key;

        // Get no. nodes for the face
        key.size = static_cast<int>(face_nodes.size());

        // Sort the nodes
        std::copy(face_nodes.cbegin(), face_nodes.cend(), key.data.begin());
        sort_4({key.data.begin(), key.data.begin() + key.size});

        if (map_.contains(key))
        {
            return map_.at(key);
        }
        else
        {
            return {};
        }
    }
    //=========================================================================
    std::vector<CellType> FaceMap::get_faces() const
    {
        std::vector<CellType> faces(map_.size());
        for (const auto &[k, v] : map_)
        {
            faces[v.second] = v.first;
        }
        return faces;
    }
    //=========================================================================
    std::pair<std::vector<CellType>,
              std::shared_ptr<graph::Connectivity>>
    extract_faces(const std::vector<Cell> &cells,
                  const graph::Connectivity &cell_to_node)
    {
        Timer timer;

        SFEM_CHECK_SIZES(cells.size(), cell_to_node.n_primary());

        // Compute the cell-to-face connectivity offsets
        std::vector<int> cell_n_faces(cells.size(), 0);
        for (std::size_t i = 0; i < cells.size(); i++)
        {
            cell_n_faces[i] = cell_num_faces(cells[i].type);
        }
        std::vector<int> cell_face_offsets(cells.size() + 1, 0);
        std::inclusive_scan(cell_n_faces.cbegin(),
                            cell_n_faces.cend(),
                            cell_face_offsets.begin() + 1);

        // Fill the cell-to-face connectivity array.
        std::vector<int> cell_face_array(cell_face_offsets.back(), 0);
        std::fill(cell_n_faces.begin(), cell_n_faces.end(), 0);
        FaceMap face_map;
        for (int i = 0; i < cell_to_node.n_primary(); i++)
        {
            auto cell = cells[i];
            auto cell_nodes = cell_to_node.links(i);

            for (int j = 0; j < cell_num_faces(cell.type); j++)
            {
                auto face_type = cell_face_type(cell.type, j);
                auto face_ordering = cell_face_ordering(cell.type, j);
                std::array<int, 4> face_nodes;
                for (int k = 0; k < cell_num_nodes(face_type); k++)
                {
                    face_nodes[k] = cell_nodes[face_ordering[k]];
                }
                cell_face_array[cell_face_offsets[i] + cell_n_faces[i]++] = face_map.insert({face_nodes.cbegin(),
                                                                                             face_nodes.cbegin() + cell_num_nodes(face_type)},
                                                                                            face_type);
            }
        }

        return {face_map.get_faces(),
                std::make_shared<graph::Connectivity>(std::move(cell_face_offsets),
                                                      std::move(cell_face_array))};
    }
    //=========================================================================
    std::shared_ptr<graph::Connectivity>
    face_to_node(const std::vector<Cell> &cells,
                 const graph::Connectivity &cell_to_face,
                 const graph::Connectivity &cell_to_node,
                 const IndexMap &cell_index_map)
    {
        SFEM_CHECK_SIZES(cells.size(), cell_to_face.n_primary());
        SFEM_CHECK_SIZES(cells.size(), cell_to_node.n_primary());
        SFEM_CHECK_SIZES(cells.size(), cell_index_map.n_local());

        // Face-to-cell connectivity
        auto face_to_cell = cell_to_face.invert();

        // Compute the face-to-node connectivity offsets
        std::vector<int> face_node_offsets(face_to_cell.n_primary() + 1, 0);
        for (int i = 0; i < face_to_cell.n_primary(); i++)
        {
            auto face_cells = face_to_cell.links(i);
            auto face_type = cell_face_type(cells[face_cells[0]].type,
                                            cell_to_face.relative_index(face_cells[0], i));
            face_node_offsets[i] = cell_num_nodes(face_type);
        }
        std::exclusive_scan(face_node_offsets.cbegin(),
                            face_node_offsets.cend(),
                            face_node_offsets.begin(), 0);

        // Fill the face-to-node connectivity array
        // The face's owner cell is the adjacent cell with the highest
        // global index
        std::vector<int> face_node_array(face_node_offsets.back(), 0);
        for (int i = 0; i < face_to_cell.n_primary(); i++)
        {
            auto face_cells = cell_index_map.local_to_global(face_to_cell.links(i));
            auto owner_cell_idx = cell_index_map.global_to_local(std::ranges::max(face_cells));
            auto owner_cell_type = cells[owner_cell_idx].type;
            auto owner_cell_nodes = cell_to_node.links(owner_cell_idx);
            auto face_rel_idx = cell_to_face.relative_index(owner_cell_idx, i);
            auto face_type = cell_face_type(owner_cell_type, face_rel_idx);
            auto face_ordering = cell_face_ordering(owner_cell_type, face_rel_idx);
            for (int j = 0; j < cell_num_nodes(face_type); j++)
            {
                face_node_array[face_node_offsets[i] + j] = owner_cell_nodes[face_ordering[j]];
            }
        }

        return std::make_shared<graph::Connectivity>(std::move(face_node_offsets),
                                                     std::move(face_node_array));
    }
}