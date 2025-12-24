#pragma once

#include <sfem/mesh/mesh.hpp>

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

        /// @brief Get the midpoint for a given cell
        std::array<real_t, 3> cell_midpoint(int cell_idx) const;

        /// @brief Get the volume for a given cell
        real_t cell_volume(int cell_idx) const;

        /// @brief Get the midpoint for a given facet
        std::array<real_t, 3> facet_midpoint(int facet_idx) const;

        /// @brief Get the area vector for a given facet
        geo::Vec3 facet_area_vec(int facet_idx) const;

        /// @brief Get the indices of the two adjacent cells at a facet
        std::array<int, 2> facet_adjacent_cells(int facet_idx) const;

        /// @brief Get the distances between the adjacent cell midpoints
        /// and the facet midpoint, for a given facet
        std::array<real_t, 2> facet_cell_distances(int facet_idx) const;

        /// @brief Get the distance between the midpoints of the
        /// two adjacent cells at a given facet
        geo::Vec3 facet_intercell_distance(int facet_idx) const;

        /// @brief Get the geometric interpolation factor for a facet
        real_t facet_interp_factor(int facet_idx) const;

        /// @brief Check whether a given facet is located on the boundary
        bool is_boundary(int facet_idx) const;

        /// @brief Decompose a facet's area vector into the orthogonal and nonorthogonal parts
        std::array<geo::Vec3, 2> decompose_area_vec(int facet_idx) const;

    private:
        /// @brief The mesh
        std::shared_ptr<const mesh::Mesh> mesh_;

        /// @brief DoF-to-DoF connectivity
        std::shared_ptr<graph::Connectivity> connectivity_;

        /// @brief DoF index map
        std::shared_ptr<IndexMap> index_map_;

        /// @brief Cell midpoints
        std::vector<std::array<real_t, 3>> cell_midpoints_;

        /// @brief Cell volumes
        std::vector<real_t> cell_volumes_;

        /// @brief Facet midpoints
        std::vector<std::array<real_t, 3>> facet_midpoints_;

        /// @brief Facet area vectors
        std::vector<geo::Vec3> facet_area_vecs_;

        /// @brief Cells adjacent to each facet
        std::vector<std::array<int, 2>> facet_adjacent_cells_;

        /// @brief Distance between the adjacent cell midpoints and the facet midpoint
        /// for each facet
        std::vector<std::array<real_t, 2>> facet_cell_distances_;

        /// @brief Distance between the adjacent cell midpoints for each facet
        std::vector<geo::Vec3> facet_intercell_distances_;

        /// @brief Geometric interpolation factor for each facet
        /// @note Defined as the ratio of the distance between the midpoint
        /// of the non-owner cell and the facet midpoint, over the
        /// distance of the two adjacent cell midpoints
        std::vector<real_t> facet_interp_factor_;
    };
}