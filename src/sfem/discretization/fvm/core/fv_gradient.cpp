#include "fv_gradient.hpp"
#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/discretization/fvm/core/fv_bc.hpp>
#include <sfem/geo/utils.hpp>
#include <sfem/mesh/utils/geo_utils.hpp>
#include <sfem/mesh/utils/loop_utils.hpp>
#include <sfem/la/native/dense_matrix.hpp>
#include <sfem/la/native/vector.hpp>

namespace sfem::fvm
{
    //=============================================================================
    void green_gauss_gradient(FVField &phi)
    {
        // Quick access
        const auto V = phi.space();
        const int dim = V->mesh()->pdim();
        const int n_comp = phi.n_comp();
        const auto &bc = phi.boundary_condition();
        auto &grad = phi.grad();

        // Zero the gradient
        grad.set_all(0.0);

        // Construct gradient
        auto facet_work = [&](const mesh::Mesh &,
                              const mesh::Region &region,
                              const mesh::Cell &,
                              int facet_idx)
        {
            // Cells adjacent to the facet
            const auto [owner, neighbour] = V->facet_adjacent_cells(facet_idx);

            // Facet area vector
            const geo::Vec3 Sf = V->facet_area_vec(facet_idx);

            // Boundary facets
            if (owner == neighbour)
            {
                for (int i = 0; i < n_comp; i++)
                {
                    const BCType bc_type = bc.region_type(region.name());

                    real_t phif = 0.0;
                    if (bc_type == BCType::dirichlet)
                    {
                        phif = bc.facet_value(facet_idx, i);
                    }
                    else if (bc_type == BCType::neumann)
                    {
                        const real_t dfP = V->facet_cell_distances(facet_idx)[0];
                        phif = phi.cell_value(owner) + bc.facet_value(facet_idx, i) * dfP;
                    }
                    else if (bc_type == BCType::robin)
                    {
                        /// @todo
                    }
                    else // Zero-Neumann
                    {
                        phif = phi.cell_value(owner, i);
                    }

                    // Add gradient contributions for the adjacent cell
                    for (int j = 0; j < dim; j++)
                    {
                        grad(owner, i * dim + j) += phif * Sf(j);
                    }
                }
            }
            // Internal facets
            else
            {
                // Compute gradient componentwise
                for (int i = 0; i < n_comp; i++)
                {
                    const real_t phif = phi.facet_value(facet_idx, i);

                    // Add gradient contributions for the adjacent cells
                    for (int j = 0; j < dim; j++)
                    {
                        grad(owner, i * dim + j) += phif * Sf(j);
                        grad(neighbour, i * dim + j) -= phif * Sf(j);
                    }
                }
            }
        };
        mesh::utils::for_all_facets(*V->mesh(), facet_work);

        // Gradient was incrementally constructed - needs assembly
        grad.assemble();

        // Scale by the (inverse) cell volume
        auto cell_work = [&](const mesh::Mesh &,
                             const mesh::Region &,
                             const mesh::Cell &,
                             int cell_idx)
        {
            const real_t vol_inv = 1.0 / V->cell_volume(cell_idx);
            for (int i = 0; i < grad.block_size(); i++)
            {
                grad(cell_idx, i) *= vol_inv;
            }
        };
        mesh::utils::for_all_cells(*V->mesh(), cell_work);

        // Ensure that ghost cell gradients are also updated
        grad.update_ghosts();
    }
    //=============================================================================
    void least_squares_gradient(FVField &phi)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const auto topo = mesh->topology();
        const int dim = mesh->pdim();
        const int n_comp = phi.n_comp();
        auto &grad = phi.grad();

        // Zero the gradient
        grad.set_all(0.0);

        // Least-squares problem LHS matrix and RHS vector
        la::DenseMatrix A(dim, dim);
        la::DenseMatrix b(dim, n_comp);

        // Construct gradient
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            // Reset least-squares LHS and RHS
            A.set_all(0.0);
            b.set_all(0.0);

            // Cell midpoint
            const auto xP = V->cell_midpoint(cell_idx);

            for (const auto neighour : topo->adjacent_entities(cell_idx, dim, dim))
            {
                // Neighbour cell midpoint
                const auto xN = V->cell_midpoint(neighour);

                // Intercell distance vector
                const geo::Vec3 dPN(xP, xN);

                // Weighting factor
                const real_t w = 1.0 / dPN.mag();

                // Add contributions for this neighbour
                for (int i = 0; i < dim; i++)
                {
                    for (int j = 0; j < dim; j++)
                    {
                        A(i, j) += w * dPN(i) * dPN(j);
                    }

                    for (int k = 0; k < n_comp; k++)
                    {
                        b(i, k) += w * dPN(i) * (phi.cell_value(neighour, k) - phi.cell_value(cell_idx, k));
                    }
                }
            }

            // Solve least squares problem and compute cell gradient
            const auto grad_cell = A.invert().first * b;
            for (int i = 0; i < dim; i++)
            {
                for (int k = 0; k < n_comp; k++)
                {
                    grad(cell_idx, k * dim + i) = grad_cell(i, k);
                }
            }
        };
        mesh::utils::for_all_cells(*mesh, work);

        // Ensure that ghost cell gradients are also updated
        grad.update_ghosts();
    }
}