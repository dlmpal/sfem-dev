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
                   const std::vector<std::vector<real_t>> &cell_values,
                   const std::vector<std::string> &cell_names,
                   const std::vector<std::vector<real_t>> &node_values,
                   const std::vector<std::string> &node_names)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Check that sizes match
        SFEM_CHECK_SIZES(cell_types.size(), cell_to_node.n_primary());
        SFEM_CHECK_SIZES(cell_to_node.n_secondary(), points.size());

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

        // Node values
        if (node_values.size() > 0)
        {
            SFEM_CHECK_SIZES(node_values.size(), node_names.size());
            file << std::format("POINT_DATA {}\n", points.size());
            for (std::size_t i = 0; i < node_values.size(); i++)
            {
                SFEM_CHECK_SIZES(points.size(), node_values[i].size());
                file << std::format("SCALARS {} float\n", node_names[i]);
                file << "LOOKUP_TABLE default\n";
                for (std::size_t j = 0; j < node_values[i].size(); j++)
                {
                    file << node_values[i][j] << "\n";
                }
            }
        }

        // Cell values
        if (cell_values.size() > 0)
        {
            SFEM_CHECK_SIZES(cell_values.size(), cell_names.size());
            file << std::format("CELL_DATA {}\n", cell_to_node.n_primary());
            for (std::size_t i = 0; i < cell_values.size(); i++)
            {
                SFEM_CHECK_SIZES(cell_to_node.n_primary(), cell_values[i].size());
                file << std::format("SCALARS {} float\n", cell_names[i]);
                file << "LOOKUP_TABLE default\n";
                for (std::size_t j = 0; j < cell_values[i].size(); j++)
                {
                    file << cell_values[i][j] << "\n";
                }
            }
        }
    }
}