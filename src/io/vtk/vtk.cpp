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
               const std::vector<std::shared_ptr<const Function>> &cell_funcs,
               const std::vector<std::shared_ptr<const Function>> &node_funcs,
               VTKFileType type)
    {
        // Check that sizes match
        SFEM_CHECK_SIZES(cell_types.size(), cell_to_node.n_primary());
        SFEM_CHECK_SIZES(cell_to_node.n_secondary(), points.size());
        for (const auto &func : cell_funcs)
        {
            SFEM_CHECK_SIZES(cell_types.size(), func->n_local());
        }
        for (const auto &func : node_funcs)
        {
            SFEM_CHECK_SIZES(points.size(), func->n_local());
        }

        if (type == VTKFileType::legacy)
        {
            legacy::write_vtk(filename.replace_extension(".vtk"),
                              cell_types, cell_to_node, points,
                              cell_funcs, node_funcs);
        }
        else if (type == VTKFileType::xml)
        {
            // Create the directory under which the individual .vtu files will be placed
            std::filesystem::create_directories(filename.parent_path() / filename.stem());

            // Create the individual .vtu files (one for each process)
            xml::write_vtu(filename.parent_path() / filename.stem() / std::format("proc_{}.vtu", mpi::rank()),
                           cell_types, cell_to_node, points,
                           cell_funcs, node_funcs);

            // Root process creates the .pvtu file
            int n_procs = mpi::n_procs();
            if (mpi::rank() == mpi::root())
            {
                std::vector<std::filesystem::path> sources;
                for (int i = 0; i < n_procs; i++)
                {
                    sources.push_back(filename.stem() / std::format("proc_{}.vtu", i));
                }

                xml::write_pvtu(filename.replace_extension(".pvtu"), sources, cell_funcs, node_funcs);
            }
        }
    }
    //=============================================================================
    void write(std::filesystem::path filename,
               const mesh::Mesh &mesh,
               const std::vector<std::shared_ptr<const Function>> &cell_funcs,
               const std::vector<std::shared_ptr<const Function>> &node_funcs,
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
              cell_funcs, node_funcs, type);
    }
    //=============================================================================
    void write(const std::filesystem::path &filename,
               const std::vector<std::shared_ptr<const fem::FEFunction>> &funcs,
               VTKFileType type)
    {
        if (funcs.empty())
        {
            SFEM_ERROR("Cannot create VTK file with 0 FEFunctions\n");
        }

        // Quick access
        const auto fe_space = funcs.front()->space();
        const auto mesh = fe_space->mesh();
        const auto topology = mesh->topology();
        int dim = topology->dim();

        // Cells types
        // Cells have order equal to that of the finite element space
        std::vector<int> cell_types(topology->n_entities(dim));
        for (int i = 0; i < topology->n_entities(dim); i++)
        {
            cell_types[i] = cell_type_to_vtk(topology->entity(i, dim).type,
                                             fe_space->order());
        }

        // Cast FEFunctions to Functions
        std::vector<std::shared_ptr<const Function>> funcs_;
        for (const auto &func : funcs)
        {
            if (func->space()->name() != fe_space->name())
            {
                SFEM_ERROR("All FEFunctions must belong to the same space\n");
            }
            funcs_.push_back(func);
        }

        write(filename, cell_types, *fe_space->connectivity()[0],
              fe_space->dof_points(), {}, funcs_, type);
    }
}