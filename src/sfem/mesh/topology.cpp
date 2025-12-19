#include "topology.hpp"
#include <sfem/mesh/partition.hpp>
#include <sfem/mesh/utils/face_utils.hpp>
#include <sfem/mesh/utils/edge_utils.hpp>
#include <sfem/base/error.hpp>

namespace sfem::mesh
{
    //=============================================================================
    Topology::Topology(std::vector<Cell> &&cells,
                       std::shared_ptr<const IndexMap> cell_index_map,
                       std::shared_ptr<const graph::Connectivity> cell_to_node)
        : cells_(std::move(cells))
    {
        SFEM_CHECK_SIZES(cells_.size(), cell_index_map->n_local());
        SFEM_CHECK_SIZES(cells_.size(), cell_to_node->n_primary());

        // Get the topological dimension, and make require that
        // all cells are of the same dimension
        dim_ = cell_dim(cells_.front().type);
        if (!std::all_of(cells_.cbegin(),
                         cells_.cend(),
                         [&](const Cell &cell)
                         {
                             return cell_dim(cell.type) == dim_;
                         }))
        {
            SFEM_ERROR(std::format("Not all cells are of the same topological dimension ({})\n", dim_));
        }

        // Cell-to-entity connectivities
        {
            // Cell-to-cell
            connectivity_[dim_][dim_] = std::make_shared<graph::Connectivity>(cell_to_node->primary_to_primary());

            // Cell-to-face
            if (dim_ > 2)
            {
                std::tie(std::ignore, connectivity_[dim_][2]) = utils::extract_faces(cells_, *cell_to_node);
            }

            // Cell-to-edge
            if (dim_ > 1)
            {
                connectivity_[dim_][1] = utils::extract_edges(cells_, *cell_to_node);
            }

            // Cell-to-node
            connectivity_[dim_][0] = cell_to_node;
        }

        // Partition all cell-to-entity connectivites
        // using the existing cell partition, while also
        // obtaining the corresponding enity index maps
        index_map_[dim_] = cell_index_map;
        for (int i = 0; i < dim_; i++)
        {
            std::tie(index_map_[i],
                     connectivity_[dim_][i]) = mesh::create_entity_partition(*index_map_[dim_],
                                                                             *connectivity_[dim_][i]);
        }

        // Create the facets
        const auto cell_to_facet = connectivity_[dim_][dim_ - 1];
        facets_.resize(cell_to_facet->n_secondary());
        for (int i = 0; i < cell_to_facet->n_primary(); i++)
        {
            auto cell_facets = cell_to_facet->links(i);
            for (int j = 0; j < cell_to_facet->n_links(i); j++)
            {
                int facet_idx = cell_facets[j];

                // Depending on the topological dimension,
                // the facets are either faces (dim=3),
                // edges (dim=2), or nodes (dim=1). For dim=3,
                // the face type (triangle or quadrilateral), must be
                // specified
                if (dim_ - 1 == 2)
                {
                    facets_[facet_idx].type = cell_face_type(cells_[i].type, j);
                }
                else if (dim_ - 1 == 1)
                {
                    facets_[facet_idx].type = CellType::line;
                }
                else
                {
                    facets_[facet_idx].type = CellType::point;
                }
            }
        }

        // Face-to-entity connectivities
        if (dim_ > 2)
        {
            // Face-to-node
            connectivity_[2][0] = utils::face_to_node(cells_,
                                                      *connectivity_[dim_][2],
                                                      *connectivity_[dim_][0],
                                                      *index_map_[dim_]);

            // Face-to-edge
            connectivity_[2][1] = utils::face_to_edge(cells_, facets_,
                                                      *connectivity_[dim_][1],
                                                      *connectivity_[dim_][0],
                                                      *connectivity_[2][0]);

            // Face-to-face
            connectivity_[2][2] = std::make_shared<graph::Connectivity>(connectivity_[2][0]->primary_to_primary());

            // Face-to-cell
            connectivity_[2][dim_] = std::make_shared<graph::Connectivity>(connectivity_[dim_][2]->invert());
        }

        // Edge-to-entity connectivities
        if (dim_ > 1)
        {
            // Edge-to-node
            connectivity_[1][0] = utils::edge_to_node(cells_,
                                                      *connectivity_[dim_][1],
                                                      *connectivity_[dim_][0],
                                                      *index_map_[dim_]);

            // Edge-to-edge
            connectivity_[1][1] = std::make_shared<graph::Connectivity>(connectivity_[1][0]->primary_to_primary());

            // Edge-to-face
            connectivity_[1][2] = std::make_shared<graph::Connectivity>(connectivity_[2][1]->invert());

            // Edge-to-cell
            connectivity_[1][dim_] = std::make_shared<graph::Connectivity>(connectivity_[dim_][1]->invert());
        }

        // Node-to-entity connectivities
        {
            // Node-to-cell
            connectivity_[0][dim_] = std::make_shared<graph::Connectivity>(connectivity_[dim_][0]->invert());

            // Node-to-node
            if (dim_ == 1)
            {
                /// @todo Clean-up
                int n_nodes = connectivity_[dim_][0]->n_secondary();
                std::vector<int> offsets(n_nodes + 1, 0);
                std::ranges::iota(offsets, 0);
                std::vector<int> array(n_nodes, 0);
                std::ranges::iota(array, 0);
                connectivity_[0][0] = std::make_shared<graph::Connectivity>(std::move(offsets),
                                                                            std::move(array));
            }
            else
            {
                connectivity_[0][0] = std::make_shared<graph::Connectivity>(connectivity_[0][dim_]->primary_to_primary());
            }

            // Node-to-edge
            if (dim_ > 1)
            {
                connectivity_[0][1] = std::make_shared<graph::Connectivity>(connectivity_[1][0]->invert());
            }

            // Node-to-face
            if (dim_ > 2)
            {
                connectivity_[0][2] = std::make_shared<graph::Connectivity>(connectivity_[2][0]->invert());
            }
        }

        // Set the facet tags
        for (int i = 0; i < n_entities(dim_ - 1); i++)
        {
            const int owner_cell_idx = entity_owner(i, dim_ - 1);
            facets_[i].tag = cells_[owner_cell_idx].tag;
        }
    }
    //=============================================================================
    const std::vector<Cell> &Topology::cells() const
    {
        return cells_;
    }
    //=============================================================================
    const std::vector<Cell> &Topology::facets() const
    {
        return facets_;
    }
    //=============================================================================
    std::shared_ptr<const graph::Connectivity> Topology::connectivity(int dim1, int dim2) const
    {
        if (dim1 < 0 or dim1 > dim_)
        {
            /// @todo error
        }
        if (dim2 < 0 or dim2 > dim_)
        {
            /// @todo error
        }
        return connectivity_[dim1][dim2];
    }
    //=============================================================================
    std::shared_ptr<const IndexMap> Topology::entity_index_map(int dim) const
    {
        if (dim < 0 or dim > dim_)
        {
            /// @todo error
        }
        return index_map_[dim];
    }
    //=============================================================================
    int Topology::dim() const
    {
        return dim_;
    }
    //=============================================================================
    int Topology::n_entities(int dim) const
    {
        return connectivity_[dim][0]->n_primary();
    }
    //=============================================================================
    Cell Topology::entity(int entity_idx, int dim) const
    {
        if (dim > dim_ or dim < 0)
        {
            /// @todo error
        }

        if (dim == dim_)
        {
            return cells_[entity_idx];
        }
        else if (dim == dim_ - 1)
        {
            return facets_[entity_idx];
        }
        else if (dim == dim_ - 2 and dim_ == 3)
        {
            return Cell{.tag = 0, .type = CellType::line};
        }
        else
        {
            return Cell{.tag = 0, .type = CellType::point};
        }
    }
    //=============================================================================
    std::span<const int> Topology::adjacent_entities(int entity_idx, int dim1, int dim2) const
    {
        if (dim1 > dim_ or dim1 < 0)
        {
            /// @todo error
        }

        if (dim2 > dim_ or dim2 < 0)
        {
            /// @todo error
        }

        return connectivity_[dim1][dim2]->links(entity_idx);
    }
    //=============================================================================
    int Topology::entity_rel_idx(int entity_idx, int dim1,
                                 int adjacent_entity_local_idx, int dim2) const
    {
        return connectivity_[dim1][dim2]->relative_index(entity_idx,
                                                         adjacent_entity_local_idx);
    }
    //=============================================================================
    int Topology::entity_owner(int entity_idx, int dim) const
    {
        if (dim == dim_)
        {
            return entity_idx;
        }
        const auto cell_im = index_map_[dim_];
        auto entity_cells = cell_im->local_to_global(connectivity_[dim][dim_]->links(entity_idx));
        return cell_im->global_to_local(std::ranges::max(entity_cells));
    }
    //=============================================================================
    std::array<int, 2> Topology::facet_adjacent_cells(int facet_local_idx) const
    {
        auto facet_cells = connectivity_[dim_ - 1][dim_]->links(facet_local_idx);
        auto owner_cell = entity_owner(facet_local_idx, dim_ - 1);
        if (facet_cells.size() == 1)
        {
            return {owner_cell, owner_cell};
        }
        else
        {
            std::array<int, 2> adjacent_cells = {facet_cells[0],
                                                 facet_cells[1]};
            if (adjacent_cells[0] != owner_cell)
            {
                std::swap(adjacent_cells[0], adjacent_cells[1]);
            }
            return adjacent_cells;
        }
    }
    //=============================================================================
    void Topology::set_cell_tag(int cell_local_idx, int tag)
    {
        SFEM_CHECK_INDEX(cell_local_idx, n_entities(dim_));
        cells_[cell_local_idx].tag = tag;
    }
    //=============================================================================
    void Topology::set_facet_tag(int facet_local_idx, int tag)
    {
        SFEM_CHECK_INDEX(facet_local_idx, n_entities(dim_ - 1));
        facets_[facet_local_idx].tag = tag;
    }
}