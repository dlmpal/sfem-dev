#pragma once

#include "topology.hpp"
#include <memory>
#include <vector>

namespace sfem::mesh
{
    /// @brief This class is able to fully describe a computational mesh,
    /// as it is composed by a topology (i.e. the connectivity relations
    /// between mesh entities) and a geometry (i.e. geometrical info about
    /// the mesh entities)
    class Mesh
    {
    public:
        /// @brief Create a Mesh
        /// @param topology The mesh topology (cells, faces etc.)
        /// @param points Node x,y,z coordinates
        /// @param regions The mesh regions
        Mesh(std::shared_ptr<Topology> topology,
             std::vector<std::array<real_t, 3>> &&points,
             std::vector<Region> &&regions);

        // Disable copying
        Mesh(const Mesh &) = delete;
        Mesh &operator=(const Mesh &) = delete;

        // Move constructor and assignment
        Mesh(Mesh &&) = default;
        Mesh &operator=(Mesh &&mesh) = default;

        /// @brief Get a reference to the mesh topology
        std::shared_ptr<Topology> topology();

        /// @brief Get a reference to the mesh topology (const version)
        std::shared_ptr<const Topology> topology() const;

        /// @brief Get a reference to the mesh node coordinates
        std::vector<std::array<real_t, 3>> &points();

        /// @brief Get a reference to the mesh node coordinates (const version)
        const std::vector<std::array<real_t, 3>> &points() const;

        /// @brief Get a reference to the mesh regions
        const std::vector<Region> &regions() const;

        /// @brief Get the physical dimension of the mesh
        /// @note Might differ from the topological dimension
        int pdim() const;

        /// @brief Topological dimension of the mesh
        int tdim() const;

        /// @brief Get the coordinates of an entity's nodes
        /// @param entity_idx The entity's index
        /// @param dim The entity's topological dimension
        /// @return The entity's points
        std::vector<std::array<real_t, 3>>
        entity_points(int entity_idx, int dim) const;

        /// @brief Get a region by its name
        Region get_region_by_name(const std::string &region_name) const;

        /// @brief Get all cells belonging to a region
        std::vector<std::pair<Cell, int>>
        region_cells(const std::string &region_name) const;

        /// @brief Get all faces belonging to a region
        std::vector<std::pair<Cell, int>>
        region_facets(const std::string &region_name) const;

    private:
        /// @brief Mesh topology
        std::shared_ptr<Topology> topology_;

        /// @brief Node coordinates
        std::vector<std::array<real_t, 3>> points_;

        /// @brief Regions
        std::vector<Region> regions_;

        /// @brief Physical dimension
        int dim_;
    };
}
