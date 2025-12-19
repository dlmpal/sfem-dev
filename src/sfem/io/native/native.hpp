#pragma once

#include <sfem/mesh/sfem_mesh.hpp>
#include <filesystem>

namespace sfem::io
{
    /// @brief Read and partition a mesh from file
    /// @param directory The mesh directory
    /// @param partition_criterion The criterion for constructing the cell-to-cell connectivity
    /// @param partitioner_type The partitioner type
    /// @return The mesh. When using multiple processes, returns part of the mesh
    /// belonging to this process
    /// @note Inside the mesh directory, the files with the following names are expected:
    /// [cells, boundary_facets, regions, points]
    /// If any file is missing, calling this function will result in error
    std::shared_ptr<mesh::Mesh>
    read_mesh(const std::filesystem::path &directory,
              mesh::PartitionCriterion partition_criterion = mesh::PartitionCriterion::shared_node,
              graph::partition::PartitionerType partitioner_type = graph::partition::PartitionerType::metis);

    /// @brief Write a mesh to file
    /// @param directory The mesh directory
    /// @param mesh The mesh
    /// @note Not MPI collective
    void write_mesh(const std::filesystem::path &directory, const mesh::Mesh &mesh);
}

namespace sfem::io::native
{
    /// @brief Write the cells and the cell-to-node connectivity to file
    /// @param filename The filename
    /// @param topology The mesh topology
    void write_cells(const std::filesystem::path &filename,
                     const mesh::Topology &topology);

    /// @brief Read the cells and the cell-to-node connectivity
    /// @param filename The filename
    /// @param cell_im The cell indexmap. Only local cells are stored
    /// @return The cells, the cell-to-node connectivity, and
    /// a map for mapping nodes from global to local indexing
    /// @note If an empty indexmap is passed, i.e. an indexmap
    /// with 0 owned indices, all cells are stored
    std::tuple<std::vector<mesh::Cell>,
               std::shared_ptr<graph::Connectivity>,
               std::unordered_map<int, int>>
    read_cells(const std::filesystem::path &filename,
               const IndexMap &cell_im);

    /// @brief Write the boundary facets to file
    /// @param filename The filename
    /// @param topology The mesh topology
    void write_boundary_facets(const std::filesystem::path &filename,
                               const mesh::Topology &topology);

    /// @brief Read the boundary facets from file
    /// @param filename The filename
    /// @param topology The mesh topology
    void read_boundary_facets(const std::filesystem::path &filename,
                              mesh::Topology &topology);

    /// @brief Write the xyz coordinates of mesh nodes to file
    /// @param filename The filename
    /// @param points The nodal coordinates
    void write_points(const std::filesystem::path &filename,
                      const std::vector<std::array<real_t, 3>> &points);

    /// @brief Read the xyz coordinates of mesh nodes from file
    /// @param filename The filename
    /// @param global_to_local The node global-to-local mapping.
    /// Only the coordinates of local nodes are stored.
    /// @return The xyz coordinates for each mesh node
    std::vector<std::array<real_t, 3>>
    read_points(const std::filesystem::path &filename,
                const std::unordered_map<int, int> &global_to_local);

    /// @brief Write the mesh regions to file
    /// @param filename The filename
    /// @param regions The regions
    void write_regions(const std::filesystem::path &filename,
                       const std::vector<mesh::Region> &regions);

    /// @brief Read the mesh regions from file
    /// @param filename The filename
    /// @return The mesh regions
    std::vector<mesh::Region>
    read_regions(const std::filesystem::path &filename);
}