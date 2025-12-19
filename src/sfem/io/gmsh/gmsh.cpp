#include "gmsh.hpp"
#include "utils.hpp"
#include <sfem/mesh/utils/face_utils.hpp>
#include <sfem/base/error.hpp>
#include <sfem/base/timer.hpp>
#include <numeric>
#include <format>
#include <fstream>

namespace sfem::io::gmsh
{
    //=============================================================================
    std::shared_ptr<mesh::Mesh> read(const std::filesystem::path &filename)
    {
        Timer timer;

        std::ifstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Lambda to skip forward in the file
        auto skip = [&file](std::size_t n_pos) -> std::size_t
        {std::string buffer; std::size_t i=0; 
        for(;i < n_pos; i++) {file >> buffer;} return i; };

        // Keep track of position in file
        std::size_t file_idx = 0;

        // Skip the first 4 lines
        file_idx += skip(6);

        // $PhysicalNames section
        // Read number of regions
        int n_regions;
        file >> n_regions;
        file_idx++;

        // Store regions into vector
        std::vector<mesh::Region> regions;
        for (int i = 0; i < n_regions; i++)
        {
            std::string name;
            int dim, tag;
            file >> dim;
            file >> tag;
            file >> name;
            std::erase(name, '"');
            regions.emplace_back(name, tag, dim);
            file_idx += 3;
        }
        int mesh_dim = std::ranges::max(regions,
                                        [](mesh::Region &lhs, mesh::Region &rhs)
                                        { return lhs.dim() < rhs.dim(); })
                           .dim();

        // Skip the following two lines
        file_idx += skip(2);

        // $Nodes section
        // Read number of nodes
        int n_nodes;
        file >> n_nodes;
        file_idx++;

        // Store nodal positions into vector
        std::vector<std::array<real_t, 3>> points(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            int node_idx;
            file >> node_idx;
            file >> points[i][0];
            file >> points[i][1];
            file >> points[i][2];
            file_idx += 4;
        }

        // Skip the following two lines
        file_idx += skip(2);

        // $Elements section
        // Read number total of elements/cells
        int n_elements;
        file >> n_elements;
        file_idx++;

        // Cells vector
        std::vector<mesh::Cell> cells;

        // Cell-to-node connectivity offsets
        std::vector<int> cell_node_offsets;

        // Cell-to-node connectivity array
        std::vector<int> cell_node_array;

        // Boundary facet map
        // Instead of using an EdgeMap for 2D meshes,
        // the FaceMap can also work, albeit less
        // efficiently
        mesh::utils::FaceMap bfacet_map;

        // Boundary facet tags
        std::vector<int> bfacet_tags;

        // Lambda to read the $Elements section of the gmsh file
        // During the first call, the function stores the number of nodes
        // per cell, and fills the cells vector.
        // During the second call, the function fills the connectivity
        // array of the cell-to-node connectivity.
        auto read_cells = [&](bool first_call)
        {
            // Number of nodes per cell
            std::vector<int> cell_n_nodes;

            // Current cell index
            int cell_idx = 0;

            // Current boundary facet index
            int bfacet_idx = 0;

            for (int ele_idx = 0; ele_idx < n_elements; ele_idx++)
            {
                int gmsh_ele_idx, gmsh_ele_type;
                file >> gmsh_ele_idx;
                file >> gmsh_ele_type;

                int n_tags;
                file >> n_tags;

                int physical_tag, elementary_tag;
                file >> physical_tag;
                file >> elementary_tag;

                auto type = gmsh::gmsh_type_to_native(gmsh_ele_type, gmsh_ele_idx);

                if (mesh::cell_dim(type) == mesh_dim)
                {
                    if (first_call)
                    {
                        cells.emplace_back(mesh::Cell{.tag = physical_tag,
                                                      .type = type});
                        cell_n_nodes.emplace_back(mesh::cell_num_nodes(type));
                        skip(mesh::cell_num_nodes(type));
                    }
                    else
                    {
                        int node_idx;
                        int offset = cell_node_offsets[cell_idx];
                        for (int i = 0; i < mesh::cell_num_nodes(cells[cell_idx].type); i++)
                        {
                            file >> node_idx;
                            cell_node_array[offset++] = node_idx - 1;
                        }
                    }
                    cell_idx++;
                }
                else
                {
                    if (first_call)
                    {
                        std::array<int, 4> facet_nodes;
                        int node_idx;
                        for (int i = 0; i < mesh::cell_num_nodes(type); i++)
                        {
                            file >> node_idx;
                            facet_nodes[i] = node_idx - 1;
                        }
                        bfacet_map.insert({facet_nodes.cbegin(),
                                           facet_nodes.cbegin() + mesh::cell_num_nodes(type)},
                                          type);
                        bfacet_tags.emplace_back(physical_tag);
                    }
                    else
                    {
                        skip(mesh::cell_num_nodes(type));
                    }
                    bfacet_idx++;
                }
            }

            if (first_call)
            {
                // Resize cells vector so that capacity=size
                cells.shrink_to_fit();

                // Compute the cell-to-node connectivity offsets
                cell_node_offsets.resize(cells.size() + 1, 0);
                std::inclusive_scan(cell_n_nodes.cbegin(),
                                    cell_n_nodes.cend(),
                                    cell_node_offsets.begin() + 1);

                // Initialize the cell-to-node connectivity array
                cell_node_array.resize(cell_node_offsets.back(), 0);

                // The file is closed and re-opened at the ($Elements) section
                file.close();
                file.open(filename);
                skip(file_idx);
            }
        };

        // First pass, store cells and compute cell-to-node offsets
        read_cells(true);

        // Second pass, fill the connectivity array
        read_cells(false);

        // Create the mesh topology
        auto topology = std::make_shared<mesh::Topology>(std::move(cells),
                                                         std::make_shared<IndexMap>(static_cast<int>(cells.size())),
                                                         std::make_shared<graph::Connectivity>(std::move(cell_node_offsets),
                                                                                               std::move(cell_node_array)));

        // Assign correct region tag to boundary facets
        int dim = topology->dim(); // Topological dimension
        for (int i = 0; i < topology->n_entities(dim - 1); i++)
        {
            auto facet = topology->entity(i, dim - 1);

            // Check if the facet is a boundary facet
            // If it is, assign the correct region tag
            if (auto adjacent_cells = topology->facet_adjacent_cells(i);
                adjacent_cells[0] == adjacent_cells[1])
            {
                auto bfacet_nodes = topology->adjacent_entities(i, dim - 1, 0);
                auto bfacet = bfacet_map.at({bfacet_nodes.cbegin(),
                                             bfacet_nodes.cbegin() + mesh::cell_num_nodes(facet.type)});
                if (bfacet.has_value())
                {
                    topology->set_facet_tag(i, bfacet_tags[bfacet.value().second]);
                }
            }
        }

        return std::make_shared<mesh::Mesh>(topology,
                                            std::move(points),
                                            std::move(regions));
    }
}