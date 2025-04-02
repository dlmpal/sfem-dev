#pragma once

#include "cell.hpp"
#include "region.hpp"
#include "../graph/connectivity.hpp"
#include "../parallel/index_map.hpp"
#include <memory>

namespace sfem::mesh
{
    /// @brief This class is tasked with storing connectivity information
    /// about all mesh entities (i.e. cells, faces, edges and nodes).
    class Topology
    {
    public:
        Topology(std::vector<Cell> &&cells,
                 std::shared_ptr<const IndexMap> cell_index_map,
                 std::shared_ptr<const graph::Connectivity> cell_to_node);

        /// @brief Get the mesh cells
        const std::vector<Cell> &cells() const;

        /// @brief Get the mesh facets
        const std::vector<Cell> &facets() const;

        /// @brief Get the connectivity between two sets of mesh entities
        std::shared_ptr<const graph::Connectivity> connectivity(int dim1, int dim2) const;

        /// @brief Get the index map for entities of a given dimension
        std::shared_ptr<const IndexMap> entity_index_map(int dim) const;

        /// @brief The topological dimension
        /// i.e. the highest dimension of the entities
        int dim() const;

        /// @brief Get the number of entities of a given dimension
        int n_entities(int dim) const;

        /// @brief Get an entity of given dimension using its index
        Cell entity(int entity_idx, int dim) const;

        /// @brief For an entity of dimension dim1, get its adjacent entities
        /// of dimension dim2
        std::span<const int>
        adjacent_entities(int entity_idx, int dim1, int dim2) const;

        /// @brief Get the relative index of an entity of dimension dim2,
        /// adjacent to one of dimension dim1
        int entity_rel_idx(int entity_idx, int dim1,
                           int adjacent_entity_idx, int dim2) const;

        /// @brief Get the owner cell for an entity of dimension dim
        /// @note The owner cell is adjacent cell with the highest global index
        /// @note If dim is the topological dimension, i.e. the entity is
        /// itself a cell, returns the cell's index
        int entity_owner(int entity_idx, int dim) const;

        /// @brief Get the cells adjacent to the facet. The first adjacent
        /// cell is always the owner cell. The orientation of the facet
        /// matches that cell's perspective. If the facet is boundary,
        /// i.e. adjacent to only one cell, then the owner cell is returned
        /// twice.
        std::array<int, 2> facet_adjacent_cells(int facet_idx) const;

        /// @brief Set the integer tag for a mesh cell
        void set_cell_tag(int cell_idx, int tag);

        /// @brief Set the integer tag for a mesh facet
        void set_facet_tag(int facet_idx, int tag);

    private:
        /// @brief Internal cells
        std::vector<Cell> cells_;

        /// @brief Internal and boundary facets
        std::vector<Cell> facets_;

        /// @brief Connectivities between all mesh entities
        std::array<std::array<std::shared_ptr<const graph::Connectivity>, 4>, 4> connectivity_;

        /// @brief Index maps for all mesh entities
        std::array<std::shared_ptr<const IndexMap>, 4> index_map_;

        /// @brief Topological dimension
        int dim_;
    };
}
