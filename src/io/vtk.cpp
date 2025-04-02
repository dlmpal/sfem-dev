#include "vtk.hpp"
#include "../base/error.hpp"
#include <map>
#include <ranges>
#include <format>
#include <fstream>
#include <filesystem>

namespace sfem::io::vtk
{
    //=============================================================================
    static int cell_type_to_vtk(mesh::CellType cell_type, int order)
    {
        static const std::map<std::pair<mesh::CellType, int>, int> to_vtk = {
            {{mesh::CellType::point, 1}, 1},

            {{mesh::CellType::line, 1}, 3},
            {{mesh::CellType::line, 2}, 21},
            {{mesh::CellType::line, 3}, 68},

            {{mesh::CellType::triangle, 1}, 5},
            {{mesh::CellType::triangle, 2}, 22},
            {{mesh::CellType::triangle, 3}, 69},

            {{mesh::CellType::quad, 1}, 9},
            {{mesh::CellType::quad, 2}, 23},
            {{mesh::CellType::quad, 3}, 70},

            {{mesh::CellType::tet, 1}, 10},
            {{mesh::CellType::tet, 2}, 24},
            {{mesh::CellType::tet, 3}, 71},

            {{mesh::CellType::hex, 1}, 12},
            {{mesh::CellType::hex, 2}, 25},
            {{mesh::CellType::hex, 3}, 72},

            {{mesh::CellType::prism, 1}, 13}

        };

        if (!to_vtk.contains({cell_type, order}))
        {
            SFEM_ERROR(std::format("Cell type {} of degree {} cannot be converted to VTK equivalent\n",
                                   static_cast<int>(cell_type), order));
        }

        return to_vtk.at({cell_type, order});
    }
    //=============================================================================
    void write(const std::filesystem::path &filename,
               const std::vector<mesh::Cell> &cells,
               const std::vector<int> &cell_orders,
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
        SFEM_CHECK_SIZES(cells.size(), cell_to_node.n_primary());
        SFEM_CHECK_SIZES(cells.size(), cell_orders.size());
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
        file << std::format("CELL_TYPES {}\n", cells.size());
        for (const auto &[cell, order] : std::ranges::zip_view(cells, cell_orders))
        {
            file << vtk::cell_type_to_vtk(cell.type, order) << "\n";
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
            file << std::format("CELL_DATA {}\n", cells.size());
            for (std::size_t i = 0; i < cell_values.size(); i++)
            {
                SFEM_CHECK_SIZES(cells.size(), cell_values[i].size());
                file << std::format("SCALARS {} float\n", cell_names[i]);
                file << "LOOKUP_TABLE default\n";
                for (std::size_t j = 0; j < cell_values[i].size(); j++)
                {
                    file << cell_values[i][j] << "\n";
                }
            }
        }
    }
    //=============================================================================
    void write(const std::filesystem::path &filename,
               const mesh::Mesh &mesh)
    {
        int dim = mesh.topology()->dim();
        std::vector<int> cell_orders(mesh.topology()->cells().size(), 1);
        write(filename,
              mesh.topology()->cells(),
              cell_orders,
              *mesh.topology()->connectivity(dim, 0),
              mesh.points());
    }
    //=============================================================================
    void write(const std::filesystem::path &filename,
               const fem::FESpace &fe_space,
               const std::vector<real_t> &values)
    {
        std::vector<int> cell_orders(fe_space.mesh()->topology()->cells().size(),
                                     fe_space.order());

        std::vector<std::string> comp_names = fe_space.components();
        int n_comp = fe_space.n_comp();
        int n_dof = static_cast<int>(values.size()) / n_comp;
        std::vector<std::vector<real_t>> comp_values(n_comp);
        for (int i = 0; i < n_comp; i++)
        {
            comp_values[i].resize(n_dof);
        }
        for (int i = 0; i < n_dof; i++)
        {
            for (int j = 0; j < n_comp; j++)
            {
                comp_values[j][i] = values[i * n_comp + j];
            }
        }

        write(filename,
              fe_space.mesh()->topology()->cells(),
              cell_orders,
              *fe_space.connectivity()[0],
              fe_space.dof_points(),
              {}, {},
              comp_values, comp_names);
    }
}