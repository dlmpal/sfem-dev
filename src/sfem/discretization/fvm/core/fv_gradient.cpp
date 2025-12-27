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
        auto &grad = phi.grad();

        // Zero the gradient
        grad.set_all(0.0);

        // Construct gradient
        auto facet_work = [&](const mesh::Mesh &,
                              const mesh::Region &,
                              const mesh::Cell &,
                              int facet_idx)
        {
            const auto [owner, neighbour] = V->facet_adjacent_cells(facet_idx);
            const geo::Vec3 Sf = V->facet_area_vec(facet_idx);

            for (int i = 0; i < n_comp; i++)
            {
                const real_t phif = phi.facet_value(facet_idx, i);
                for (int j = 0; j < dim; j++)
                {
                    grad(owner, i * dim + j) += phif * Sf(j);
                    if (owner != neighbour)
                    {
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
                        int owner)
        {
            // Reset least-squares LHS and RHS
            A.set_all(0.0);
            b.set_all(0.0);

            // Cell midpoint
            const auto xP = V->cell_midpoint(owner);

            for (const auto facet : topo->adjacent_entities(owner, dim, dim - 1))
            {
                const auto adjacent_cells = V->facet_adjacent_cells(facet);
                const auto neighbour = owner == adjacent_cells[0] ? adjacent_cells[1]
                                                                  : adjacent_cells[0];

                // Owner-to-neighbour distance vector
                geo::Vec3 dPN;
                if (owner == neighbour)
                {
                    dPN = 2.0 * geo::Vec3(xP, V->facet_midpoint(facet));
                }
                else
                {
                    dPN = geo::Vec3(xP, V->cell_midpoint(neighbour));
                }

                // Weighting factor
                const real_t w = 1.0 / dPN.mag();

                // Add contributions for this neighbour
                for (int i = 0; i < dim; i++)
                {
                    for (int j = 0; j < dim; j++)
                    {
                        A(i, j) += w * dPN(i) * dPN(j);
                    }

                    if (owner == neighbour)
                    {
                        for (int k = 0; k < n_comp; k++)
                        {
                            b(i, k) += 2.0 * w * dPN(i) * (phi.facet_value(facet, k) - phi.cell_value(owner, k));
                        }
                    }
                    else
                    {
                        for (int k = 0; k < n_comp; k++)
                        {
                            b(i, k) += w * dPN(i) * (phi.cell_value(neighbour, k) - phi.cell_value(owner, k));
                        }
                    }
                }
            }

            // Solve least squares problem and compute cell gradient
            const auto grad_cell = A.invert().first * b;
            for (int i = 0; i < dim; i++)
            {
                for (int k = 0; k < n_comp; k++)
                {
                    grad(owner, k * dim + i) = grad_cell(i, k);
                }
            }
        };
        mesh::utils::for_all_cells(*mesh, work);

        // Ensure that ghost cell gradients are also updated
        grad.update_ghosts();
    }
}