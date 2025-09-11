#pragma once

#include "../../mesh/mesh.hpp"

namespace sfem::fvm
{
    /// @brief A FVSpace stores and readily provides
    /// geometric information often used in finite volume
    /// discretizations, e.g. facet area and normal vector,
    /// along with DoF information, similar to FESpace
    /// for finite element problems
    class FVSpace
    {
    public:
        /// @brief Create a FVSpace for a given mesh
        /// @param mesh The mesh
        FVSpace(std::shared_ptr<const mesh::Mesh> mesh);

        /// @brief Get the mesh
        std::shared_ptr<const mesh::Mesh> mesh() const;

        /// @brief Get the DoF-to-DoF connectivity
        std::shared_ptr<const graph::Connectivity> connectivity() const;

        /// @brief Get the DoF index map
        std::shared_ptr<const IndexMap> index_map() const;

        /// @brief Get the volume for a given cell
        real_t cell_volume(int cell_idx) const;

        /// @brief Get the area for a given facet
        real_t facet_area(int facet_idx) const;

        /// @brief Get the midpoint for a given cell
        std::array<real_t, 3> cell_midpoint(int cell_idx) const;

        /// @brief Get the midpoint for a given facet
        std::array<real_t, 3> facet_midpoint(int facet_idx) const;

        /// @brief Get the normal vector for a given facet
        /// @note Normalized
        geo::Vec3 facet_normal(int facet_idx) const;

        /// @brief Get the distance between the midpoints of the
        /// two adjacent cells at a given facet
        geo::Vec3 intercell_distance(int facet_idx) const;

        /// @brief Get the distances between the adjacent cell midpoints
        /// and the facet midpoint, for a given facet
        std::array<real_t, 2> facet_cell_distances(int facet_idx) const;

        /// @brief Get the indices of the two adjacent cells at a facet
        std::array<int, 2> facet_adjacent_cells(int facet_idx) const;

        /// @brief Check whether a given facet is located on the boundary
        bool is_boundary(int facet_idx) const;

        /// @brief Compute a distance weighted average at the facet of two
        /// adjacent cells, given the value of each cell
        real_t compute_facet_value(int facet_idx,
                                   real_t value1,
                                   real_t value2,
                                   bool harmonic = false) const;

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