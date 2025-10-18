#include "legacy.hpp"
#include "utils.hpp"
#include "../../base/error.hpp"
#include <map>
#include <format>
#include <fstream>
#include <filesystem>
#include <ranges>

namespace sfem::io::vtk::legacy
{
    //=============================================================================
    void write_vtk(const std::filesystem::path &filename,
                   const std::vector<int> &cell_types,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<std::shared_ptr<const Field>> &cell_fields,
                   const std::vector<std::shared_ptr<const Field>> &node_fields)
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
            std::vector<int> cell_nodes;
            cell_nodes.assign(cell_to_node.links(i).cbegin(),
                              cell_to_node.links(i).cend());
            cell_node_ordering_to_vtk(cell_types[i], cell_nodes);
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
        for (const auto &func : cell_fields)
        {
            for (const auto &comp : func->components())
            {
                file << std::format("SCALARS {} float\n", comp);
                file << "LOOKUP_TABLE default\n";
                const int comp_idx = func->comp_idx(comp);
                for (int i = 0; i < func->n_local(); i++)
                {
                    file << (*func)(i, comp_idx) << "\n";
                }
            }
        }

        // Point data
        file << std::format("POINT_DATA {}\n", cell_to_node.n_secondary());
        for (const auto &func : node_fields)
        {
            for (const auto &comp : func->components())
            {
                file << std::format("SCALARS {} float\n", comp);
                file << "LOOKUP_TABLE default\n";
                const int comp_idx = func->comp_idx(comp);
                for (int i = 0; i < func->n_local(); i++)
                {
                    file << (*func)(i, comp_idx) << "\n";
                }
            }
        }
    }
}