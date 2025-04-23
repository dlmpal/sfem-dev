#include "vtk.hpp"
#include "legacy.hpp"
#include "xml.hpp"
#include "../../parallel/mpi.hpp"
#include <format>

namespace sfem::io::vtk
{
    //=============================================================================
    void write(std::filesystem::path filename,
               const std::vector<mesh::Cell> &cells,
               const std::vector<int> &cell_orders,
               const graph::Connectivity &cell_to_node,
               const std::vector<std::array<real_t, 3>> &points,
               const std::vector<std::vector<real_t>> &cell_values,
               const std::vector<std::string> &cell_names,
               const std::vector<std::vector<real_t>> &node_values,
               const std::vector<std::string> &node_names,
               VTKFileType type)
    {
        if (type == VTKFileType::legacy)
        {
            legacy::write_vtk(filename.replace_extension(".vtk"),
                              cells, cell_orders,
                              cell_to_node, points,
                              cell_values, cell_names,
                              node_values, node_names);
        }
        else if (type == VTKFileType::xml)
        {
            // Create the directory under which the individual .vtu files will be placed
            std::filesystem::create_directories(filename.parent_path() / filename.stem());

            // Create the individual .vtu files (one for each process)
            xml::write_vtu(filename.parent_path() / filename.stem() / std::format("proc_{}.vtu", mpi::rank()),
                           cells, cell_orders,
                           cell_to_node, points,
                           cell_values, cell_names,
                           node_values, node_names);

            // Root process creates the .pvtu file
            int n_procs = mpi::n_procs();
            if (mpi::rank() == mpi::root())
            {
                std::vector<std::filesystem::path> sources;
                for (int i = 0; i < n_procs; i++)
                {
                    sources.push_back(filename.stem() / std::format("proc_{}.vtu", i));
                }

                xml::write_pvtu(filename.replace_extension(".pvtu"), sources, cell_names, node_names);
            }
        }
    }
    //=============================================================================
    void write(std::filesystem::path filename,
               const mesh::Mesh &mesh,
               VTKFileType type)
    {
        // Quick access
        const auto topology = mesh.topology();
        int dim = topology->dim();

        // All cells are of order 1 (no FESpace involved)
        std::vector<int> cell_orders(mesh.topology()->cells().size(), 1);

        write(filename,
              topology->cells(), cell_orders,
              *topology->connectivity(dim, 0), mesh.points(),
              {}, {}, {}, {},
              type);
    }
    //=============================================================================
    void write(const std::filesystem::path &filename,
               const fem::FESpace &fe_space,
               const std::vector<real_t> &values,
               VTKFileType type)
    {
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

        // All cells have the order of the space
        std::vector<int> cell_orders(fe_space.mesh()->topology()->cells().size(),
                                     fe_space.order());

        write(filename,
              fe_space.mesh()->topology()->cells(), cell_orders,
              *fe_space.connectivity()[0], fe_space.dof_points(),
              {}, {},
              comp_values, fe_space.components(),
              type);
    }
}