#pragma once

#include "../../mesh/mesh.hpp"

namespace sfem::fvm
{
    /// @brief
    class FVSpace
    {
    public:
        FVSpace(std::shared_ptr<const mesh::Mesh> mesh);

        std::shared_ptr<const mesh::Mesh> mesh() const;

        std::shared_ptr<const graph::Connectivity> connectivity() const;

        std::shared_ptr<const IndexMap> index_map() const;

        real_t cell_volume(int cell_idx) const;

        real_t facet_area(int facet_idx) const;

        std::array<real_t, 3> cell_midpoint(int cell_idx) const;

        std::array<real_t, 3> facet_midpoint(int facet_idx) const;

        geo::Vec3 facet_normal(int facet_idx) const;

        geo::Vec3 intercell_distance(int facet_idx) const;

        std::array<real_t, 2> facet_cell_distances(int facet_idx) const;

        std::array<int, 2> facet_adjacent_cells(int facet_idx) const;

        bool is_boundary(int facet_idx) const;

    private:
        /// @brief The mesh
        std::shared_ptr<const mesh::Mesh> mesh_;

        /// @brief DoF-to-DoF connectivity
        std::shared_ptr<graph::Connectivity> connectivity_;

        /// @brief DoF index map
        std::shared_ptr<IndexMap> index_map_;

        /// @brief Cell volumes
        std::vector<real_t> cell_volumes_;

        /// @brief Facet areas
        std::vector<real_t> facet_areas_;

        /// @brief Cell centers
        std::vector<std::array<real_t, 3>> cell_midpoints_;

        /// @brief Facet centers
        std::vector<std::array<real_t, 3>> facet_midpoints_;

        /// @brief Facet normal vectors
        std::vector<geo::Vec3> facet_normals_;

        /// @brief Distance between adjacent cell centers at each facet
        std::vector<geo::Vec3> intercell_distances_;

        /// @brief Distance between adjacent cell centers and the facet center
        std::vector<std::array<real_t, 2>> facet_cell_distances_;

        /// @brief Cells adjacent to each facet
        std::vector<std::array<int, 2>> facet_adjacent_cells_;
    };
}