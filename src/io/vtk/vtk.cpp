#include "vtk.hpp"
#include "utils.hpp"
#include "legacy.hpp"
#include "xml.hpp"
#include "../../parallel/mpi.hpp"
#include <format>

namespace sfem::io::vtk
{
    //=============================================================================
    void write(std::filesystem::path filename,
               const std::vector<int> &cell_types,
               const graph::Connectivity &cell_to_node,
               const std::vector<std::array<real_t, 3>> &points,
               const std::vector<std::pair<std::string, std::span<real_t>>> &cell_data,
               const std::vector<std::pair<std::string, std::span<real_t>>> &node_data,
               VTKFileType type)
    {
        // Check that sizes match
        SFEM_CHECK_SIZES(cell_types.size(), cell_to_node.n_primary());
        SFEM_CHECK_SIZES(cell_to_node.n_secondary(), points.size());
        for (const auto &[name, data] : cell_data)
        {
            SFEM_CHECK_SIZES(cell_types.size(), data.size());
        }
        for (const auto &[name, data] : node_data)
        {
            SFEM_CHECK_SIZES(points.size(), data.size());
        }

        if (type == VTKFileType::legacy)
        {
            legacy::write_vtk(filename.replace_extension(".vtk"),
                              cell_types, cell_to_node, points,
                              cell_data, node_data);
        }
        else if (type == VTKFileType::xml)
        {
            // Create the directory under which the individual .vtu files will be placed
            std::filesystem::create_directories(filename.parent_path() / filename.stem());

            // Create the individual .vtu files (one for each process)
            xml::write_vtu(filename.parent_path() / filename.stem() / std::format("proc_{}.vtu", mpi::rank()),
                           cell_types, cell_to_node, points,
                           cell_data, node_data);

            // Root process creates the .pvtu file
            int n_procs = mpi::n_procs();
            if (mpi::rank() == mpi::root())
            {
                std::vector<std::filesystem::path> sources;
                for (int i = 0; i < n_procs; i++)
                {
                    sources.push_back(filename.stem() / std::format("proc_{}.vtu", i));
                }

                xml::write_pvtu(filename.replace_extension(".pvtu"), sources, cell_data, node_data);
            }
        }
    }
    //=============================================================================
    void write(std::filesystem::path filename,
               const mesh::Mesh &mesh,
               const std::vector<std::pair<std::string, std::span<real_t>>> &cell_data,
               const std::vector<std::pair<std::string, std::span<real_t>>> &node_data,
               VTKFileType type)
    {
        // Quick access
        const auto topology = mesh.topology();
        int dim = topology->dim();
        const auto &cell_to_node = *topology->connectivity(dim, 0);

        // All cells are of order 1 (no finite element space involved)
        std::vector<int> cell_types(topology->n_entities(dim));
        for (int i = 0; i < topology->n_entities(dim); i++)
        {
            cell_types[i] = cell_type_to_vtk(topology->entity(i, dim).type, 1);
        }

        write(filename, cell_types, cell_to_node, mesh.points(),
              cell_data, node_data, type);
    }
    //=============================================================================
    void write(const std::filesystem::path &filename,
               const fem::FESpace &fe_space,
               const std::vector<real_t> &values,
               VTKFileType type)
    {
        /// @todo Try to avoid the copies
        // Split components into separate vectors
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
        std::vector<std::pair<std::string, std::span<real_t>>> node_data;
        for (int i = 0; i < n_comp; i++)
        {
            node_data.push_back({fe_space.components()[i], comp_values[i]});
        }

        // Quick access
        const auto mesh = fe_space.mesh();
        const auto topology = mesh->topology();
        int dim = topology->dim();

        // Cells types
        // Cells have order equal to that of the finite element space
        std::vector<int> cell_types(topology->n_entities(dim));
        for (int i = 0; i < topology->n_entities(dim); i++)
        {
            cell_types[i] = cell_type_to_vtk(topology->entity(i, dim).type,
                                             fe_space.order());
        }

        write(filename, cell_types, *fe_space.connectivity()[0],
              fe_space.dof_points(), {}, node_data, type);
    }
}