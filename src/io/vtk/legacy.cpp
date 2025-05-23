#include "legacy.hpp"
#include "../../base/error.hpp"
#include <map>
#include <format>
#include <fstream>
#include <filesystem>

namespace sfem::io::vtk::legacy
{
    //=============================================================================
    void write_vtk(const std::filesystem::path &filename,
                   const std::vector<int> &cell_types,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<std::pair<std::string, std::span<real_t>>> &cell_data,
                   const std::vector<std::pair<std::string, std::span<real_t>>> &node_data)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // File header
        file << "# vtk DataFile Version 2.0\n";
        file << "SFEM\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";

        // Points
        file << std::format("POINTS {} float\n", points.size());
        for (std::size_t i = 0; i < points.size(); i++)
        {
            file << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
        }

        // Cells
        file << std::format("CELLS {} {}\n",
                            cell_to_node.n_primary(),
                            cell_to_node.n_primary() + cell_to_node.n_links());
        for (int i = 0; i < cell_to_node.n_primary(); i++)
        {
            file << cell_to_node.n_links(i) << " ";
            auto cell_nodes = cell_to_node.links(i);
            for (int j = 0; j < cell_to_node.n_links(i); j++)
            {
                file << cell_nodes[j] << " ";
            }
            file << "\n";
        }

        // Cell types
        file << std::format("CELL_TYPES {}\n", cell_to_node.n_primary());
        for (int i = 0; i < cell_to_node.n_primary(); i++)
        {
            file << cell_types[i] << "\n";
        }

        // Cell data
        file << std::format("CELL_DATA {}\n", cell_to_node.n_primary());
        for (const auto &[name, data] : cell_data)
        {
            file << std::format("SCALARS {} float\n", name);
            file << "LOOKUP_TABLE default\n";
            for (std::size_t i = 0; i < data.size(); i++)
            {
                file << data[i] << "\n";
            }
        }

        // Point data
        file << std::format("POINT_DATA {}\n", cell_to_node.n_secondary());
        for (const auto &[name, data] : node_data)
        {
            file << std::format("SCALARS {} float\n", name);
            file << "LOOKUP_TABLE default\n";
            for (std::size_t i = 0; i < data.size(); i++)
            {
                file << data[i] << "\n";
            }
        }
    }
}