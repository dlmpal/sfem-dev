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
               const std::vector<IOField> &cell_fields,
               const std::vector<IOField> &node_fields,
               VTKFileType type)
    {
        // Check that sizes match
        SFEM_CHECK_SIZES(cell_types.size(), cell_to_node.n_primary());
        SFEM_CHECK_SIZES(cell_to_node.n_secondary(), points.size());

        if (type == VTKFileType::legacy)
        {
            legacy::write_vtk(filename.replace_extension(".vtk"),
                              cell_types, cell_to_node, points,
                              cell_fields, node_fields);
        }
        else if (type == VTKFileType::xml)
        {
            // Create the directory under which the individual .vtu files will be placed
            std::filesystem::create_directories(filename.parent_path() / filename.stem());

            // Create the individual .vtu files (one for each process)
            xml::write_vtu(filename.parent_path() / filename.stem() / std::format("proc_{}.vtu", mpi::rank()),
                           cell_types, cell_to_node, points,
                           cell_fields, node_fields);

            // Root process creates the .pvtu file
            int n_procs = mpi::n_procs();
            if (mpi::rank() == mpi::root())
            {
                std::vector<std::filesystem::path> sources;
                for (int i = 0; i < n_procs; i++)
                {
                    sources.push_back(filename.stem() / std::format("proc_{}.vtu", i));
                }

                xml::write_pvtu(filename.replace_extension(".pvtu"), sources, cell_fields, node_fields);
            }
        }
    }
    //=============================================================================
    void write(std::filesystem::path filename,
               const mesh::Mesh &mesh,
               const std::vector<IOField> &cell_fields,
               const std::vector<IOField> &node_fields,
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
              cell_fields, node_fields, type);
    }
    //=============================================================================
    void write(const std::filesystem::path &filename,
               const std::vector<fvm::FVField> &fields,
               VTKFileType type)
    {
        std::vector<IOField> cell_fields;
        for (const auto &field : fields)
        {
            cell_fields.emplace_back(field.components(), [field](int idx, int comp_idx)
                                     { return field.cell_value(idx, comp_idx); });
        }

        write(filename, *fields.front().space()->mesh(), cell_fields, {}, type);
    }
    //=============================================================================
    void write(const std::filesystem::path &filename,
               const std::vector<fem::FEField> &node_fields,
               const std::vector<fem::FEField> &cell_fields,
               VTKFileType type)
    {
        if (node_fields.empty() && cell_fields.empty())
        {
            SFEM_ERROR("Cannot create VTK file with 0 FEFields\n");
        }

        // Quick access
        const auto fe_space = node_fields.front().space();
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

        std::vector<IOField> node_fields_;
        for (const auto &field : node_fields)
        {
            if (field.space()->name() != fe_space->name())
            {
                SFEM_ERROR("All FEFunctions must belong to the same space\n");
            }
            node_fields_.emplace_back(field.components(), [field](int idx, int comp_idx)
                                      { return field.dof_values()(idx, comp_idx); });
        }

        std::vector<IOField> cell_fields_;
        for (const auto &field : cell_fields)
        {
            cell_fields_.emplace_back(field.components(), [field](int idx, int comp_idx)
                                      { return field.dof_values()(idx, comp_idx); });
        }

        if (fe_space->order() == 0)
        {
            write(filename, cell_types, *topology->connectivity(dim, 0),
                  mesh->points(), node_fields_, {}, type);
        }
        else
        {
            write(filename, cell_types, *fe_space->connectivity()[0],
                  fe_space->dof_points(), cell_fields_, node_fields_, type);
        }
    }
}