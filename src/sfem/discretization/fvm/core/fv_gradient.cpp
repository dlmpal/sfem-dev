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
    void green_gauss_gradient(FVField phi)
    {
        // Quick access
        const FVSpace &V = *phi.space();
        const int dim = V.mesh()->pdim();
        const int n_comp = phi.n_comp();
        const FVBC &bc = phi.boundary_condition();
        la::Vector &grad = phi.grad();

        // Zero the gradient
        grad.set_all(0.0);

        // Construct gradient
        auto facet_work = [&](const mesh::Mesh &,
                              const mesh::Region &region,
                              const mesh::Cell &,
                              int facet_idx)
        {
            // Cells adjacent to the facet
            const auto [owner, neighbour] = V.facet_adjacent_cells(facet_idx);

            // Facet area and normal vector
            const real_t Af = V.facet_area(facet_idx);
            const auto nf = V.facet_normal(facet_idx);

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
                        const real_t d = V.facet_cell_distances(facet_idx)[0];
                        phif = phi.cell_value(owner) + bc.facet_value(facet_idx, i) * d;
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
                        grad(owner, i * dim + j) += phif * Af * nf(j);
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
                        grad(owner, i * dim + j) += phif * Af * nf(j);
                        grad(neighbour, i * dim + j) -= phif * Af * nf(j);
                    }
                }
            }
        };
        mesh::utils::for_all_facets(*V.mesh(), facet_work);

        // Gradient was incrementally constructed - needs assembly
        grad.assemble();

        // Scale by the (inverse) cell volume
        auto cell_work = [&](const mesh::Mesh &,
                             const mesh::Region &,
                             const mesh::Cell &,
                             int cell_idx)
        {
            const real_t vol_inv = 1.0 / V.cell_volume(cell_idx);
            for (int i = 0; i < grad.block_size(); i++)
            {
                grad(cell_idx, i) *= vol_inv;
            }
        };
        mesh::utils::for_all_cells(*V.mesh(), cell_work);

        // Ensure that ghost cell gradients are also updated
        grad.update_ghosts();
    }
    //=============================================================================
    // void least_squares_gradient(const FVField &phi, FVField &grad)
    // {
    //     // Zero the gradient
    //     grad.set_all(0.0);

    //     // Finite volume space, mesh and topology
    //     const auto V = phi.space();
    //     const auto mesh = V->mesh();
    //     const auto topo = mesh->topology();

    //     // Mesh dimension and number of components
    //     const int dim = mesh->pdim();
    //     const int n_comp = phi.n_comp();

    //     // Least-squares problem LHS matrix and RHS vector
    //     la::DenseMatrix A(dim, dim);
    //     la::DenseMatrix b(dim, n_comp);

    //     // Construct gradient
    //     auto work = [&](const mesh::Mesh &,
    //                     const mesh::Region &,
    //                     const mesh::Cell &,
    //                     int cell_idx)
    //     {
    //         // Reset least-squares LHS and RHS
    //         A.set_all(0.0);
    //         b.set_all(0.0);

    //         // Cell midpoint
    //         const auto r1 = V->cell_midpoint(cell_idx);

    //         for (const auto neighour : topo->adjacent_entities(cell_idx, dim, dim))
    //         {
    //             // Neighbour cell midpoint
    //             const auto r2 = V->cell_midpoint(neighour);

    //             // Intercell distance vector
    //             const geo::Vec3 d12(r1, r2);

    //             // Weighting factor
    //             const real_t w = 1.0 / d12.mag();

    //             // Add contributions for this neighbour
    //             for (int i = 0; i < dim; i++)
    //             {
    //                 for (int j = 0; j < dim; j++)
    //                 {
    //                     A(i, j) += w * d12(i) * d12(j);
    //                 }

    //                 for (int k = 0; k < n_comp; k++)
    //                 {
    //                     b(i, k) += w * d12(i) * (phi(neighour, k) - phi(cell_idx, k));
    //                 }
    //             }
    //         }

    //         // Solve least squares problem and compute cell gradient
    //         const auto grad_cell = A.invert().first * b;
    //         for (int i = 0; i < dim; i++)
    //         {
    //             for (int k = 0; k < n_comp; k++)
    //             {
    //                 grad(cell_idx, k * dim + i) = grad_cell(i, k);
    //             }
    //         }
    //     };
    //     mesh::utils::for_all_cells(*mesh, work);

    //     // Ensure that ghost cell gradients are also updated
    //     grad.update_ghosts();
    // }
}