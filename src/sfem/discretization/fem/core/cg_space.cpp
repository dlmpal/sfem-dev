#include "cg_space.hpp"
#include <sfem/discretization/fem/core/elements/fe_factory.hpp>
#include <sfem/mesh/partition.hpp>
#include <format>
#include <numeric>

namespace sfem::fem
{
    //=============================================================================
    CGSpace::CGSpace(std::shared_ptr<const mesh::Mesh> mesh, int order)
        : FESpace(mesh, order, std::format("CG({})", order))
    {
        // Quick access to mesh topology
        auto topology = mesh->topology();
        int dim = topology->dim();

        // Compute the cell-to-dof connectivity offsets
        std::vector<int> cell_dof_offsets(topology->n_entities(dim) + 1, 0);
        for (int i = 0; i < topology->n_entities(dim); i++)
        {
            cell_dof_offsets[i] = dof::cell_num_dof(topology->entity(i, dim).type,
                                                    order_);
        }
        std::exclusive_scan(cell_dof_offsets.cbegin(),
                            cell_dof_offsets.cend(),
                            cell_dof_offsets.begin(), 0);

        // Store the first DoF belonging to each subentity (nodes, edges, faces)
        std::array<std::vector<int>, 3> subentity_dof;
        for (int i = 0; i < dim; i++)
        {
            subentity_dof[i].resize(topology->n_entities(i), -1);
        }

        // Fill the cell-to-dof connectivity array
        std::vector<int> cell_dof_array(cell_dof_offsets.back(), 0);
        int n_dof = 0;
        for (int i = 0; i < topology->n_entities(dim); i++)
        {
            int offset = cell_dof_offsets[i];
            auto cell_type = topology->entity(i, dim).type;
            auto cell_nodes = topology->adjacent_entities(i, dim, 0);

            // Corner-node DoF
            if (order_ > 0)
            {
                for (int node_idx : cell_nodes)
                {
                    if (subentity_dof[0][node_idx] < 0)
                    {
                        subentity_dof[0][node_idx] = n_dof++;
                    }
                    cell_dof_array[offset++] = subentity_dof[0][node_idx];
                }
            }

            // Edge DoF
            if (topology->dim() > 1 && order_ > 0)
            {
                int n_dof_edge = dof::cell_num_internal_dof(mesh::CellType::line, order_);
                auto cell_edges = topology->adjacent_entities(i, dim, 1);
                for (std::size_t j = 0; j < cell_edges.size(); j++)
                {
                    if (subentity_dof[1][cell_edges[j]] < 0)
                    {
                        subentity_dof[1][cell_edges[j]] = n_dof;
                        n_dof += n_dof_edge;
                    }

                    for (int k = 0; k < n_dof_edge; k++)
                    {
                        cell_dof_array[offset++] = subentity_dof[1][cell_edges[j]] + k;
                    }

                    // Enforce same orientation across different cells
                    if (i != topology->entity_owner(cell_edges[j], 1))
                    {
                        std::reverse(cell_dof_array.begin() + offset - n_dof_edge,
                                     cell_dof_array.begin() + offset);
                    }
                }
            }

            // Face DoF
            if (topology->dim() > 2 && order_ > 0)
            {
                /// @todo Correctly enumerate face DoF in 3D
                /// Enforce same orientation by cell global idx
            }

            // Volume DoF
            for (int j = 0; j < dof::cell_num_internal_dof(cell_type, order_); j++)
            {
                cell_dof_array[offset++] = n_dof++;
            }
        }

        // Create the DoF partition (and the cell-to-DoF connectivity) using the cell partition
        std::tie(index_map_, connectivity_[0]) = mesh::create_entity_partition(*topology->entity_index_map(dim),
                                                                               graph::Connectivity(std::move(cell_dof_offsets),
                                                                                                   std::move(cell_dof_array)));

        // Compute the DoF-to-DoF connectivity
        connectivity_[1] = std::make_shared<graph::Connectivity>(connectivity_[0]->invert().primary_to_primary(1, true));

        // Create the element collection
        for (int cell_type = 0; cell_type < static_cast<int>(mesh::CellType::n_cell_types); cell_type++)
        {
            fe_collection_[cell_type] = create_nodal_element(static_cast<mesh::CellType>(cell_type), order_);
        }
    }
}