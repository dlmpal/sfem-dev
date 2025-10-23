// Create a uniform 2D or 3D Cartesian mesh

#include "sfem.hpp"
#include "argparse.hpp"
#include <numeric>

using namespace sfem;
using namespace mesh;

int main(int argc, char *argv[])
{
    initialize(argc, argv, "cart-mesh");

    argparse::ArgParser parser;
    parser.add_argument(argparse::Argument("d", true));
    parser.add_argument(argparse::Argument("Nx", true));
    parser.add_argument(argparse::Argument("Ny", true));
    parser.add_argument(argparse::Argument("Nz", false).value<int>(1));
    parser.add_argument(argparse::Argument("x-low", true));
    parser.add_argument(argparse::Argument("y-low", true));
    parser.add_argument(argparse::Argument("z-low", false).value<real_t>(0.0));
    parser.add_argument(argparse::Argument("x-high", true));
    parser.add_argument(argparse::Argument("y-high", true));
    parser.add_argument(argparse::Argument("z-high", false).value<real_t>(1.0));
    parser.add_argument(argparse::Argument("dir", false).value("mesh"));
    parser.parse_args(argc, argv);

    const int dim = parser.get_argument("d")->value<int>();
    if (dim < 2 || dim > 3)
    {
        SFEM_ERROR(std::format("Invalid dimension: {}\n", dim));
    }

    // Enforce that 2D meshes are "flat"
    if (dim == 2)
    {
        parser.get_argument("Nz")->value(1);
        parser.get_argument("z-low")->value(0.0);
        parser.get_argument("z-high")->value(1.0);
    }

    // Number of cells per direction
    /// @todo Check that these are valid
    const int Nx = parser.get_argument("Nx")->value<int>();
    const int Ny = parser.get_argument("Ny")->value<int>();
    const int Nz = parser.get_argument("Nz")->value<int>();

    // Number of nodes per direction
    const int nx = Nx + 1;
    const int ny = Ny + 1;
    const int nz = dim == 3 ? Nz + 1 : 1;

    // Total number of cells and nodes
    const int n_cells = Nx * Ny * Nz;
    const int n_nodes = nx * ny * nz;

    // Domain length per direction
    /// @todo Check that these are valid (!=0)
    const real_t x_low = parser.get_argument("x-low")->value<real_t>();
    const real_t y_low = parser.get_argument("y-low")->value<real_t>();
    const real_t z_low = parser.get_argument("z-low")->value<real_t>();
    const real_t x_high = parser.get_argument("x-high")->value<real_t>();
    const real_t y_high = parser.get_argument("y-high")->value<real_t>();
    const real_t z_high = parser.get_argument("z-high")->value<real_t>();

    // Cell size per direction
    const real_t dx = (x_high - x_low) / static_cast<real_t>(Nx);
    const real_t dy = (y_high - y_low) / static_cast<real_t>(Ny);
    const real_t dz = (z_high - z_low) / static_cast<real_t>(Nz);

    // Cell type and number of nodes per cell
    const CellType cell_type = dim == 3 ? CellType::hexahedron : CellType::quadrilateral;
    const int n_nodes_cell = cell_num_nodes(cell_type);

    // Create the cell array
    const int internal_tag = 1;
    std::vector<Cell> cells(n_cells, Cell{.tag = internal_tag, .type = cell_type});

    // Create cell-to-node offsets
    std::vector<int> cell_node_offsets(n_cells + 1, n_nodes_cell);
    std::exclusive_scan(cell_node_offsets.cbegin(),
                        cell_node_offsets.cend(),
                        cell_node_offsets.begin(), 0);

    // Create the cell-to-node connectivity array
    std::vector<int> cell_node_array(cell_node_offsets.back());
    for (int k = 0; k < Nz; k++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                std::array<int, 8> cell_nodes;
                if (dim == 3)
                {
                    cell_nodes[0] = k * ny * nx + j * nx + i;
                    cell_nodes[1] = k * ny * nx + j * nx + i + 1;
                    cell_nodes[2] = k * ny * nx + (j + 1) * nx + i + 1;
                    cell_nodes[3] = k * ny * nx + (j + 1) * nx + i;
                    cell_nodes[4] = (k + 1) * ny * nx + j * nx + i;
                    cell_nodes[5] = (k + 1) * ny * nx + j * nx + i + 1;
                    cell_nodes[6] = (k + 1) * ny * nx + (j + 1) * nx + i + 1;
                    cell_nodes[7] = (k + 1) * ny * nx + (j + 1) * nx + i;
                }
                else
                {
                    cell_nodes[0] = j * nx + i;
                    cell_nodes[1] = j * nx + i + 1;
                    cell_nodes[2] = (j + 1) * nx + i + 1;
                    cell_nodes[3] = (j + 1) * nx + i;
                }
                const int cell_idx = k * Ny * Nx + j * Nx + i;
                std::copy(cell_nodes.cbegin(),
                          cell_nodes.cbegin() + n_nodes_cell,
                          cell_node_array.begin() + cell_node_offsets[cell_idx]);
            }
        }
    }

    // Create the cell-to-node connectivity
    auto cell_to_node = std::make_shared<graph::Connectivity>(std::move(cell_node_offsets),
                                                              std::move(cell_node_array));

    // Create the topology
    auto topology = std::make_shared<Topology>(std::move(cells),
                                               std::make_shared<IndexMap>(n_cells),
                                               cell_to_node);

    const int left = 0;
    const int right = 1;
    const int top = 2;
    const int bottom = 3;
    const int front = 4;
    const int back = 5;

    const std::array<int, 6> boundary_tags =
        {
            2, // Left
            3, // Right
            4, // Top
            5, // Bottom
            6, // Front
            7  // Back
        };

    std::array<int, 6> boundary_rel_idx;
    if (dim == 3)
    {
        boundary_rel_idx[left] = 0;
        boundary_rel_idx[right] = 5;
        boundary_rel_idx[bottom] = 4;
        boundary_rel_idx[top] = 1;
        boundary_rel_idx[front] = 2;
        boundary_rel_idx[back] = 3;
    }
    else
    {
        boundary_rel_idx[left] = 3;
        boundary_rel_idx[right] = 1;
        boundary_rel_idx[bottom] = 0;
        boundary_rel_idx[top] = 2;
    }

    // Calculate the number of boundary facets
    const int n_bfacets_2d = 2 * Nx + 2 * Ny;
    const int n_bfacets_3d = 2 * Ny * Nz + 2 * Nx * Nz + 2 * Nx * Ny;
    const int n_bfacets = dim == 3 ? n_bfacets_3d : n_bfacets_2d;

    // Store boundary facet information
    std::vector<int> bfacet_owners(n_bfacets);
    std::vector<int> bfacet_rel_idx(n_bfacets);
    std::vector<int> bfacet_tag(n_bfacets);
    int bfacet_idx = 0;

    // Left boundary
    for (int k = 0; k < Nz; k++)
    {
        for (int j = 0; j < Ny; j++)
        {
            bfacet_owners[bfacet_idx] = k * Ny * Nx + j * Nx + 0;
            bfacet_rel_idx[bfacet_idx] = boundary_rel_idx[left];
            bfacet_tag[bfacet_idx++] = boundary_tags[left];
        }
    }

    // Right boundary
    for (int k = 0; k < Nz; k++)
    {
        for (int j = 0; j < Ny; j++)
        {
            bfacet_owners[bfacet_idx] = k * Ny * Nx + j * Nx + Nx - 1;
            bfacet_rel_idx[bfacet_idx] = boundary_rel_idx[right];
            bfacet_tag[bfacet_idx++] = boundary_tags[right];
        }
    }

    // Bottom boundary
    for (int k = 0; k < Nz; k++)
    {
        for (int i = 0; i < Nx; i++)
        {
            bfacet_owners[bfacet_idx] = k * Ny * Nx + 0 * Nx + i;
            bfacet_rel_idx[bfacet_idx] = boundary_rel_idx[bottom];
            bfacet_tag[bfacet_idx++] = boundary_tags[bottom];
        }
    }

    // Top boundary
    for (int k = 0; k < Nz; k++)
    {
        for (int i = 0; i < Nx; i++)
        {
            bfacet_owners[bfacet_idx] = k * Ny * Nx + (Ny - 1) * Nx + i;
            bfacet_rel_idx[bfacet_idx] = boundary_rel_idx[top];
            bfacet_tag[bfacet_idx++] = boundary_tags[top];
        }
    }

    if (dim == 3)
    {
        // Front boundary
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                bfacet_owners[bfacet_idx] = 0 * Ny * Nx + (Ny - 1) * Nx + i;
                bfacet_rel_idx[bfacet_idx] = boundary_rel_idx[front];
                bfacet_tag[bfacet_idx++] = boundary_tags[front];
            }
        }

        // Back boundary
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                bfacet_owners[bfacet_idx] = (Nz - 1) * Ny * Nx + (Ny - 1) * Nx + i;
                bfacet_rel_idx[bfacet_idx] = boundary_rel_idx[back];
                bfacet_tag[bfacet_idx++] = boundary_tags[back];
            }
        }
    }

    SFEM_CHECK_SIZES(bfacet_idx, n_bfacets);

    // Set the correct boundary facet tags
    for (int i = 0; i < n_bfacets; i++)
    {
        const int owner = bfacet_owners[i];
        const int rel_idx = bfacet_rel_idx[i];
        const int tag = bfacet_tag[i];
        const int bfacet_idx = topology->adjacent_entities(owner, dim, dim - 1)[rel_idx];
        topology->set_facet_tag(bfacet_idx, tag);
    }

    // Create the mesh points
    std::vector<std::array<real_t, 3>> points(n_nodes);
    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                const int node_idx = k * ny * nx + j * nx + i;
                points[node_idx][0] = i * dx + x_low;
                points[node_idx][1] = j * dy + y_low;
                points[node_idx][2] = k * dz + z_low;
            }
        }
    }

    // Create the mesh regions
    std::vector<Region> regions =
        {
            Region("Internal", internal_tag, dim),
            Region("Left", boundary_tags[left], dim - 1),
            Region("Right", boundary_tags[right], dim - 1),
            Region("Bottom", boundary_tags[bottom], dim - 1),
            Region("Top", boundary_tags[top], dim - 1),

        };
    if (dim == 3)
    {
        regions.emplace_back("Front", boundary_tags[front], dim - 1);
        regions.emplace_back("Back", boundary_tags[back], dim - 1);
    }

    // Create the mesh
    Mesh mesh(topology, std::move(points), std::move(regions));

    // Finally, write the mesh to file
    io::write_mesh(parser.get_argument("dir")->value<std::string>(), mesh);

    return 0;
}