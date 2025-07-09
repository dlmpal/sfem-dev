#include "assembly.hpp"
#include "../../../mesh/utils/geo_utils.hpp"

namespace sfem::fem
{
    //=============================================================================
    void assemble_matrix_cells(const FESpace &phi, const std::string &region,
                               FECellKernel kernel, MatSet mat)
    {
        // Quick access
        const auto mesh = phi.mesh();
        const int dim = mesh->physical_dim();

        // Loop over locally owned cells
        for (const auto &[cell, cell_idx] : mesh->region_cells(region))
        {
            // Skip ghost cells
            if (mesh->topology()->entity_index_map(dim)->is_ghost(cell_idx))
            {
                continue;
            }

            // Cell info
            auto cell_dof = phi.cell_dof(cell_idx);
            auto cell_points = phi.cell_dof_points(cell_idx);

            // Integration
            auto element = phi.element(cell.type);
            int n_element_dof = static_cast<int>(cell_dof.size()) * phi.n_comp();
            la::DenseMatrix cell_mat(n_element_dof, n_element_dof);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim, element->integration_rule()->point(nqpt),
                                               cell_points);
                cell_mat += kernel(cell_idx, data) * data.detJ;
            }
            mat(cell_dof, cell_dof, cell_mat.data());
        }
    }
    //=============================================================================
    void assemble_matrix_facets(const FESpace &phi, const std::string &region,
                                FEFacetKernel kernel, MatSet mat)
    {
        // Quick access
        const auto mesh = phi.mesh();
        const int dim = mesh->physical_dim();

        // Loop over locally owned facets
        for (const auto &[facet, facet_idx] : mesh->region_facets(region))
        {
            // Skip ghost facets
            if (mesh->topology()->entity_index_map(dim - 1)->is_ghost(facet_idx))
            {
                continue;
            }

            // Facet info
            auto facet_points = phi.facet_dof_points(facet_idx);
            auto facet_dof = phi.facet_dof(facet_idx);
            auto facet_normal = mesh::facet_normal(facet.type, facet_points).normalize();

            // Integration
            auto element = phi.element(facet.type);
            int n_element_dof = static_cast<int>(facet_dof.size()) * phi.n_comp();
            la::DenseMatrix facet_mat(n_element_dof, n_element_dof);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim,
                                               element->integration_rule()->point(nqpt),
                                               facet_points);

                facet_mat += kernel(facet_idx, data, facet_normal) * data.detJ;
            }
            mat(facet_dof, facet_dof, facet_mat.data());
        }
    }
    //=============================================================================
    void assemble_vec_cells(const FESpace &phi, const std::string &region,
                            FECellKernel kernel, VecSet vec)
    {
        // Quick access
        const auto mesh = phi.mesh();
        const int dim = mesh->physical_dim();

        // Loop over locally owned cells
        for (const auto &[cell, cell_idx] : mesh->region_cells(region))
        {
            // Skip ghost cells
            if (mesh->topology()->entity_index_map(dim)->is_ghost(cell_idx))
            {
                continue;
            }

            // Cell info
            auto cell_dof = phi.cell_dof(cell_idx);
            auto cell_points = phi.cell_dof_points(cell_idx);

            // Integration
            auto element = phi.element(cell.type);
            int n_element_dof = static_cast<int>(cell_dof.size()) * phi.n_comp();
            la::DenseMatrix cell_vec(n_element_dof, 1);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim,
                                               element->integration_rule()->point(nqpt),
                                               cell_points);

                cell_vec += kernel(cell_idx, data) * data.detJ;
            }
            vec(cell_dof, cell_vec.data());
        }
    }
    //=============================================================================
    void assemble_vec_facets(const FESpace &phi, const std::string &region,
                             FEFacetKernel kernel, VecSet vec)
    {
        // Quick access
        const auto mesh = phi.mesh();
        const int dim = mesh->physical_dim();

        // Loop over locally owned facets
        for (const auto &[facet, facet_idx] : mesh->region_facets(region))
        {
            // Skip ghost facets
            if (mesh->topology()->entity_index_map(dim - 1)->is_ghost(facet_idx))
            {
                continue;
            }

            // Facet info
            auto facet_points = phi.facet_dof_points(facet_idx);
            auto facet_normal = mesh::facet_normal(facet.type, facet_points).normalize();
            auto facet_dof = phi.facet_dof(facet_idx);

            // Integration
            auto element = phi.element(facet.type);
            int n_element_dof = static_cast<int>(facet_dof.size()) * phi.n_comp();
            la::DenseMatrix facet_vec(n_element_dof, 1);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim,
                                               element->integration_rule()->point(nqpt),
                                               facet_points);

                facet_vec += kernel(facet_idx, data, facet_normal) * data.detJ;
            }
            vec(facet_dof, facet_vec.data());
        }
    }
}