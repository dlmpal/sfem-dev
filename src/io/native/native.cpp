#include "native.hpp"
#include "../../parallel/mpi.hpp"
#include "../../base/timer.hpp"
#include "../../base/error.hpp"
#include <numeric>
#include <format>
#include <fstream>

namespace sfem::io
{
    //=============================================================================
    std::shared_ptr<mesh::Mesh>
    read_mesh(const std::filesystem::path &directory,
              mesh::PartitionCriterion partition_criterion,
              graph::partition::PartitionerType partitioner_type)
    {
        Timer timer;

        // Quick access to the filenames for each mesh file
        const std::unordered_map<std::string,
                                 std::filesystem::path>
            filenames = {{"cells", directory / "cells"},
                         {"boundary_facets", directory / "boundary_facets"},
                         {"regions", directory / "regions"},
                         {"points", directory / "points"}};

        // Cells, cell-to-node connectivity
        std::vector<mesh::Cell> cells;
        auto cell_to_node = std::make_shared<graph::Connectivity>();

        // Root process reads the "whole mesh",
        // i.e. all cells and their respective nodes
        if (mpi::rank() == mpi::root())
        {
            std::tie(cells,
                     cell_to_node,
                     std::ignore) = io::native::read_cells(directory / "cells",
                                                           IndexMap(0));
        }

        // Root process partitions the cells, and distributes the partition
        auto cell_im = mesh::create_cell_partition(cells,
                                                   *cell_to_node,
                                                   partition_criterion,
                                                   partitioner_type);

        // Each process now reads its respective part of the mesh
        std::unordered_map<int, int> node_global_to_local;
        std::tie(cells,
                 cell_to_node,
                 node_global_to_local) = io::native::read_cells(filenames.at("cells"),
                                                                *cell_im);

        // Keep a copy of the cell-to-node connectivity,
        // needed later to re-order the nodal coordinates
        graph::Connectivity cell_to_node_old(cell_to_node->offsets(),
                                             cell_to_node->array());

        // Create the topology
        auto topology = std::make_shared<mesh::Topology>(std::move(cells),
                                                         cell_im,
                                                         cell_to_node);

        // Assign correct region tags to boundary facets
        io::native::read_boundary_facets(filenames.at("boundary_facets"), *topology);

        // Read the regions
        auto regions = io::native::read_regions(filenames.at("regions"));

        // Store the xyz coordinates of each mesh node
        auto points_old = io::native::read_points(filenames.at("points"),
                                                  node_global_to_local);

        // Re-order the points
        int dim = topology->dim();
        std::vector<std::array<sfem::real_t, 3>> points(points_old.size());
        for (int i = 0; i < topology->n_entities(dim); i++)
        {
            auto cell_nodes = topology->adjacent_entities(i, dim, 0);
            auto cell_nodes_old = cell_to_node_old.links(i);
            for (int j = 0; j < cell_to_node_old.n_links(i); j++)
            {

                points[cell_nodes[j]] = points_old[cell_nodes_old[j]];
            }
        }

        return std::make_shared<mesh::Mesh>(topology,
                                            std::move(points),
                                            std::move(regions));
    }
    //=============================================================================
    void write_mesh(const std::filesystem::path &directory, const mesh::Mesh &mesh)
    {
        Timer timer;

        // Create the leaf directory and any parents
        std::filesystem::create_directories(directory);

        // Create the cells, bfaces, regions and points files
        native::write_cells(directory / "cells",
                            *mesh.topology());

        native::write_boundary_facets(directory / "boundary_facets",
                                      *mesh.topology());

        native::write_regions(directory / "regions",
                              mesh.regions());

        native::write_points(directory / "points",
                             mesh.points());
    }
}
namespace sfem::io::native
{
    //=============================================================================
    void write_cells(const std::filesystem::path &filename,
                     const mesh::Topology &topology)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Number of cells and nodes
        int dim = topology.dim();
        file << std::format("{} {}\n",
                            topology.n_entities(dim),
                            topology.n_entities(0));

        // Cells
        for (int i = 0; i < topology.n_entities(dim); i++)
        {
            auto cell = topology.entity(i, dim);

            // Cell region tag and type"
            file << std::format("{} {} ",
                                cell.tag,
                                static_cast<int>(cell.type));

            // Cell nodes
            for (int node_idx : topology.adjacent_entities(i, dim, 0))
            {
                file << node_idx << " ";
            }
            file << "\n";
        }
    }
    //=============================================================================
    std::tuple<std::vector<mesh::Cell>,
               std::shared_ptr<graph::Connectivity>,
               std::unordered_map<int, int>>
    read_cells(const std::filesystem::path &filename, const IndexMap &cell_im)
    {
        std::ifstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Lambda to skip forward in the file
        auto skip = [&file](std::size_t n_pos)
        {int dummy; std::size_t i=0; 
            for(;i < n_pos; i++) {file >> dummy;} return i; };

        // Read no. cells and nodes
        int n_cells;
        int n_nodes;
        file >> n_cells;
        file >> n_nodes;

        // Read cells and cell-to-node connectivity
        if (cell_im.n_owned() == 0)
        {
            // First pass, read cells and compute the
            // cell-to-node connectivity offsets
            std::vector<mesh::Cell> cells(n_cells);
            std::vector<int> cell_n_nodes(n_cells);
            for (int i = 0; i < n_cells; i++)
            {
                int tag, type;
                file >> tag;
                file >> type;

                cells[i].tag = tag;
                cells[i].type = static_cast<mesh::CellType>(type);
                cell_n_nodes[i] = mesh::cell_num_nodes(cells[i].type);

                skip(cell_n_nodes[i]);
            }
            std::vector<int> cell_node_offsets(n_cells + 1, 0);
            std::inclusive_scan(cell_n_nodes.cbegin(),
                                cell_n_nodes.cend(),
                                cell_node_offsets.begin() + 1);

            file.close();
            file.open(filename);
            skip(2);

            // Second pass, fill the cell-to-node connectivity array
            std::vector<int> cell_node_array(cell_node_offsets.back(), 0);
            for (int i = 0; i < n_cells; i++)
            {
                skip(2);
                int offset = cell_node_offsets[i];
                for (int j = 0; j < mesh::cell_num_nodes(cells[i].type); j++)
                {
                    file >> cell_node_array[offset++];
                }
            }
            auto cell_to_node = std::make_shared<graph::Connectivity>(std::move(cell_node_offsets),
                                                                      std::move(cell_node_array));

            // Create the trivial global-to-local mapping for the nodes
            std::unordered_map<int, int> node_global_to_local;
            for (int i = 0; i < n_nodes; i++)
            {
                node_global_to_local.insert({i, i});
            }

            return {cells, cell_to_node, node_global_to_local};
        }
        else
        {
            SFEM_CHECK_SIZES(n_cells, cell_im.n_global());

            // First pass, read cells and compute the
            // cell-to-node connectivity offsets
            std::vector<mesh::Cell> cells(cell_im.n_local());
            std::vector<int> cell_n_nodes(cell_im.n_local());
            for (int i = 0; i < n_cells; i++)
            {
                int tag, type, n_nodes;
                file >> tag;
                file >> type;
                n_nodes = mesh::cell_num_nodes(static_cast<mesh::CellType>(type));

                if (int cell_local_idx = cell_im.global_to_local(i); cell_local_idx >= 0)
                {
                    cells[cell_local_idx].tag = tag;
                    cells[cell_local_idx].type = static_cast<mesh::CellType>(type);
                    cell_n_nodes[cell_local_idx] = n_nodes;
                }

                skip(n_nodes);
            }
            std::vector<int> cell_node_offsets(cell_im.n_local() + 1, 0);
            std::inclusive_scan(cell_n_nodes.cbegin(),
                                cell_n_nodes.cend(),
                                cell_node_offsets.begin() + 1);

            // Close and re-open file
            file.close();
            file.open(filename);
            skip(2);

            // Second pass, fill the cell-to-node connectivity array
            std::vector<int> cell_node_array(cell_node_offsets.back(), 0);
            for (int i = 0; i < n_cells; i++)
            {
                int tag, type, n_nodes;
                file >> tag;
                file >> type;
                n_nodes = mesh::cell_num_nodes(static_cast<mesh::CellType>(type));

                if (int cell_local_idx = cell_im.global_to_local(i); cell_local_idx >= 0)
                {
                    int offset = cell_node_offsets[cell_local_idx];
                    for (int j = 0; j < mesh::cell_num_nodes(cells[cell_local_idx].type); j++)
                    {
                        file >> cell_node_array[offset++];
                    }
                }
                else
                {
                    skip(n_nodes);
                }
            }

            // Map the nodes locally to the range:
            // [0, n_nodes_local), where n_nodes_local is the number
            // of unique nodes for the cells local to this process
            std::unordered_map<int, int> global_to_local;
            int node_local_idx = 0;
            for (std::size_t i = 0; i < cell_node_array.size(); i++)
            {
                if (!global_to_local.contains(cell_node_array[i]))
                {
                    global_to_local.insert({cell_node_array[i], node_local_idx});
                    cell_node_array[i] = node_local_idx++;
                }
                else
                {
                    cell_node_array[i] = global_to_local.at(cell_node_array[i]);
                }
            }

            // Finally, create the connectivity
            auto cell_to_node = std::make_shared<graph::Connectivity>(std::move(cell_node_offsets),
                                                                      std::move(cell_node_array));

            return {cells, cell_to_node, global_to_local};
        }
    }
    //=============================================================================
    void write_boundary_facets(const std::filesystem::path &filename,
                               const mesh::Topology &topology)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Count no. boundary facets
        int dim = topology.dim();
        int n_bfacets = 0;
        for (int i = 0; i < topology.n_entities(dim - 1); i++)
        {
            if (auto adjacent_cells = topology.facet_adjacent_cells(i);
                adjacent_cells[0] == adjacent_cells[1])
            {
                n_bfacets++;
            }
        }

        // Write no. boundary facets
        file << n_bfacets << "\n";

        // Write boundary faces
        for (int i = 0; i < topology.n_entities(dim - 1); i++)
        {
            auto facet = topology.entity(i, dim - 1);
            if (auto adjacent_cells = topology.facet_adjacent_cells(i);
                adjacent_cells[0] == adjacent_cells[1])
            {
                auto owner_cell_idx = topology.entity_owner(i, dim - 1);
                auto facet_rel_idx = topology.entity_rel_idx(owner_cell_idx, dim, i, dim - 1);
                file << std::format("{} {} {}\n", facet.tag, owner_cell_idx, facet_rel_idx);
            }
        }
    }
    //=============================================================================
    void read_boundary_facets(const std::filesystem::path &filename,
                              mesh::Topology &topology)
    {
        std::ifstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Quick access
        int dim = topology.dim();
        const auto &cell_im = topology.entity_index_map(dim);
        const auto &cell_to_facet = topology.connectivity(dim, dim - 1);

        int n_bfaces;
        file >> n_bfaces;

        for (int i = 0; i < n_bfaces; i++)
        {
            int tag, cell_global_idx, facet_idx_in_cell;
            file >> tag;
            file >> cell_global_idx;
            file >> facet_idx_in_cell;

            // Assign boundary tag only if the face is locally owned
            if (int cell_local_idx = cell_im->global_to_local(cell_global_idx); cell_local_idx >= 0)
            {
                topology.set_facet_tag(cell_to_facet->links(cell_local_idx)[facet_idx_in_cell], tag);
            }
        }
    }
    //=============================================================================
    void write_points(const std::filesystem::path &filename,
                      const std::vector<std::array<real_t, 3>> &points)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Number of points
        file << points.size() << "\n";

        for (std::size_t i = 0; i < points.size(); i++)
        {
            file << points[i][0] << " ";
            file << points[i][1] << " ";
            file << points[i][2] << "\n";
        }
    }
    //=============================================================================
    std::vector<std::array<real_t, 3>>
    read_points(const std::filesystem::path &filename,
                const std::unordered_map<int, int> &global_to_local)
    {
        std::ifstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Read number of points
        int n_points;
        file >> n_points;

        // Read points
        std::vector<std::array<real_t, 3>> points(global_to_local.size());
        std::array<real_t, 3> point;
        for (int i = 0; i < n_points; i++)
        {
            file >> point[0];
            file >> point[1];
            file >> point[2];

            if (global_to_local.contains(i))
            {
                int node_local_idx = global_to_local.at(i);
                points[node_local_idx][0] = point[0];
                points[node_local_idx][1] = point[1];
                points[node_local_idx][2] = point[2];
            }
        }

        return points;
    }
    //=============================================================================
    void write_regions(const std::filesystem::path &filename, const std::vector<mesh::Region> &regions)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Number of regions
        file << regions.size() << "\n";

        for (const auto &region : regions)
        {
            // Region name, tag, dimension
            file << std::format("{} {} {}\n",
                                region.name(),
                                region.tag(),
                                region.dim());
        }
    }
    //=============================================================================
    std::vector<mesh::Region> read_regions(const std::filesystem::path &filename)
    {
        std::ifstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Read the number of regions
        int n_regions;
        file >> n_regions;

        // Read regions
        std::vector<mesh::Region> regions;
        for (auto i = 0; i < n_regions; i++)
        {
            std::string name;
            int dim, tag;
            file >> name;
            file >> tag;
            file >> dim;
            regions.emplace_back(name, tag, dim);
        }

        return regions;
    }
}