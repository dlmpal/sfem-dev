#include "assembly.hpp"
#include "../../../mesh/utils/geo_utils.hpp"
#include "../../../parallel/mpi.hpp"

namespace sfem::fem
{
    //=============================================================================
    void assemble_matrix_cells(const FEField &phi, const std::string &region,
                               FECellKernel kernel, la::MatSet mat)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();
        const int n_comp = phi.n_comp();

        // Loop over locally owned cells
        for (const auto &[cell, cell_idx] : mesh->region_cells(region))
        {
            // Skip ghost cells
            if (mesh->topology()->entity_index_map(dim)->is_ghost(cell_idx))
            {
                continue;
            }

            // Cell info
            auto cell_dof = V->cell_dof(cell_idx);
            auto cell_points = V->cell_dof_points(cell_idx);

            // Integration
            auto element = V->element(cell.type);
            int n_element_dof = static_cast<int>(cell_dof.size()) * n_comp;
            la::DenseMatrix cell_mat(n_element_dof, n_element_dof);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim, element->integration_rule()->point(nqpt),
                                               cell_points);
                cell_mat += kernel(cell_idx, data) * data.detJ;
            }
            mat(cell_dof, cell_dof, cell_mat.values());
        }
    }
    //=============================================================================
    void assemble_matrix_facets(const FEField &phi, const std::string &region,
                                FEFacetKernel kernel, la::MatSet mat)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();
        const int n_comp = phi.n_comp();

        // Loop over locally owned facets
        for (const auto &[facet, facet_idx] : mesh->region_facets(region))
        {
            // Skip ghost facets
            if (mesh->topology()->entity_index_map(dim - 1)->is_ghost(facet_idx))
            {
                continue;
            }

            // Facet info
            auto facet_points = V->facet_dof_points(facet_idx);
            auto facet_dof = V->facet_dof(facet_idx);
            auto facet_normal = mesh::facet_normal(facet.type, facet_points).normalize();

            // Integration
            auto element = V->element(facet.type);
            int n_element_dof = static_cast<int>(facet_dof.size()) * n_comp;
            la::DenseMatrix facet_mat(n_element_dof, n_element_dof);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim,
                                               element->integration_rule()->point(nqpt),
                                               facet_points);

                facet_mat += kernel(facet_idx, data, facet_normal) * data.detJ;
            }
            mat(facet_dof, facet_dof, facet_mat.values());
        }
    }
    //=============================================================================
    void assemble_vec_cells(const FEField &phi, const std::string &region,
                            FECellKernel kernel, la::VecSet vec)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();
        const int n_comp = phi.n_comp();

        // Loop over locally owned cells
        for (const auto &[cell, cell_idx] : mesh->region_cells(region))
        {
            // Skip ghost cells
            if (mesh->topology()->entity_index_map(dim)->is_ghost(cell_idx))
            {
                continue;
            }

            // Cell info
            auto cell_dof = V->cell_dof(cell_idx);
            auto cell_points = V->cell_dof_points(cell_idx);

            // Integration
            auto element = V->element(cell.type);
            int n_element_dof = static_cast<int>(cell_dof.size()) * n_comp;
            la::DenseMatrix cell_vec(n_element_dof, 1);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim,
                                               element->integration_rule()->point(nqpt),
                                               cell_points);

                cell_vec += kernel(cell_idx, data) * data.detJ;
            }
            vec(cell_dof, cell_vec.values());
        }
    }
    //=============================================================================
    void assemble_vec_facets(const FEField &phi, const std::string &region,
                             FEFacetKernel kernel, la::VecSet vec)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();
        const int n_comp = phi.n_comp();

        // Loop over locally owned facets
        for (const auto &[facet, facet_idx] : mesh->region_facets(region))
        {
            // Skip ghost facets
            if (mesh->topology()->entity_index_map(dim - 1)->is_ghost(facet_idx))
            {
                continue;
            }

            // Facet info
            auto facet_points = V->facet_dof_points(facet_idx);
            auto facet_normal = mesh::facet_normal(facet.type, facet_points).normalize();
            auto facet_dof = V->facet_dof(facet_idx);

            // Integration
            auto element = V->element(facet.type);
            int n_element_dof = static_cast<int>(facet_dof.size()) * n_comp;
            la::DenseMatrix facet_vec(n_element_dof, 1);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim,
                                               element->integration_rule()->point(nqpt),
                                               facet_points);

                facet_vec += kernel(facet_idx, data, facet_normal) * data.detJ;
            }
            vec(facet_dof, facet_vec.values());
        }
    }
    //=============================================================================
    la::DenseMatrix integrate_cells(const FESpace &phi,
                                    const std::string &region,
                                    FECellKernel kernel)
    {
        // Quick access
        const auto mesh = phi.mesh();
        const int dim = mesh->pdim();

        // Loop over locally owned cells and accumulate result
        la::DenseMatrix result(1, 1); ///< Resize if needed
        bool first_cell = true;
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
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim, element->integration_rule()->point(nqpt),
                                               cell_points);
                auto cell_result = kernel(cell_idx, data) * data.detJ;

                // It expected that the kernel does not change
                // its output shape between cells. Thus, the reshape,
                // if needed, is only done in the first cell
                if (first_cell and (cell_result.n_rows() != result.n_rows() or
                                    cell_result.n_cols() != result.n_cols()))
                {
                    result = la::DenseMatrix(cell_result.n_rows(),
                                             cell_result.n_cols());
                }
                first_cell = false;

                result += cell_result;
            }
        }

        // Sum values across processes
        for (int i = 0; i < result.n_rows(); i++)
        {
            for (int j = 0; j < result.n_cols(); j++)
            {
                result(i, j) = mpi::reduce(result(i, j), mpi::ReduceOperation::sum);
            }
        }

        return result;
    }
    //=============================================================================
    la::DenseMatrix integrate_facets(const FESpace &phi,
                                     const std::string &region,
                                     FEFacetKernel kernel)
    {
        // Quick access
        const auto mesh = phi.mesh();
        const int dim = mesh->pdim();

        // Loop over locally owned facets and accumulate result
        la::DenseMatrix result(1, 1); ///< Resize if needed
        bool first_facet = true;
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
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim,
                                               element->integration_rule()->point(nqpt),
                                               facet_points);

                auto facet_result = kernel(facet_idx, data, facet_normal) * data.detJ;

                // It expected that the kernel does not change
                // its output shape between facets. Thus, the reshape,
                // if needed, is only done in the first facet
                if (first_facet and (facet_result.n_rows() != result.n_rows() or
                                     facet_result.n_cols() != result.n_cols()))
                {
                    result = la::DenseMatrix(facet_result.n_rows(),
                                             facet_result.n_cols());
                }
                first_facet = false;

                result += facet_result;
            }
        }

        // Sum values across processes
        for (int i = 0; i < result.n_rows(); i++)
        {
            for (int j = 0; j < result.n_cols(); j++)
            {
                result(i, j) = mpi::reduce(result(i, j), mpi::ReduceOperation::sum);
            }
        }

        return result;
    }
}