#pragma once

#include "../mesh.hpp"

namespace sfem::mesh::utils
{
    template <typename F>
    concept MeshLoopFunc = std::invocable<F, const Mesh &, const Region &, const Cell &, int>;

    /// @brief Loop over the cells of a mesh
    /// @param mesh Mesh
    /// @param func Function to be executed for every cell
    /// @param skip_ghost Whether to skip ghost cells
    inline void for_all_cells(const Mesh &mesh, MeshLoopFunc auto &&func, bool skip_ghost = true)
    {
        // Dimension of cells
        const int cell_dim = mesh.pdim();

        // Loop over all regions
        for (const auto &region : mesh.regions())
        {
            // Skip boundary regions
            if (region.dim() < cell_dim)
            {
                continue;
            }

            // Loop over all cells in region
            for (const auto &[cell, cell_idx] : mesh.region_cells(region.name()))
            {
                // Skip ghost cells if required
                if (skip_ghost and mesh.topology()->entity_index_map(cell_dim)->is_ghost(cell_idx))
                {
                    continue;
                }

                // Do work
                func(mesh, region, cell, cell_idx);
            }
        }
    }

    /// @brief Loop over the facets of a mesh
    /// @param mesh Mesh
    /// @param func Function to be executed for every facet
    /// @param skip_ghost Whether to skip ghost facets
    /// @param skip_boundary Whether to skip boundary regions
    inline void for_all_facets(const Mesh &mesh, MeshLoopFunc auto &&func, bool skip_ghost = true, bool skip_boundary = false)
    {
        // Dimension of facets
        const int facet_dim = mesh.pdim() - 1;

        // Loop over all regions
        for (const auto &region : mesh.regions())
        {
            // Skip boundary regions if required
            if (skip_boundary and region.dim() < mesh.pdim())
            {
                continue;
            }

            // Loop over all facets in region
            for (const auto &[facet, facet_idx] : mesh.region_facets(region.name()))
            {
                // Skip ghost facets if required
                if (skip_ghost and mesh.topology()->entity_index_map(facet_dim)->is_ghost(facet_idx))
                {
                    continue;
                }

                // Do work
                func(mesh, region, facet, facet_idx);
            }
        }
    }
}