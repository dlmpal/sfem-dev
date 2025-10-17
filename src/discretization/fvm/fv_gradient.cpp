#include "fv_gradient.hpp"
#include "../../geo/utils.hpp"
#include "../../mesh/utils/geo_utils.hpp"
#include "../../la/native/dense_matrix.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void green_gauss_gradient(const FVFunction &phi,
                              const FVBC &bc,
                              FVFunction &grad)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();

        // Zero the gradient
        // grad.set_all(0.0, true);
        grad.set_all(0.0);

        // Loop over all facets
        for (const auto &region : mesh->regions())
        {
            for (const auto &[facet, facet_idx] : mesh->region_facets(region.name()))
            {
                // Skip ghost facets
                if (mesh->topology()->entity_index_map(mesh->pdim() - 1)->is_ghost(facet_idx))
                {
                    continue;
                }

                // Cells adjacent to the facet
                const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
                const auto &[cell_idx1, cell_idx2] = adjacent_cells;

                // Facet area and normal vector
                const real_t area = V->facet_area(facet_idx);
                const auto normal = V->facet_normal(facet_idx);

                // Boundary facets
                if (cell_idx1 == cell_idx2)
                {
                    // Boundary facet value
                    // Corresponds to zero Neumann BC by default
                    real_t phi_facet = phi(cell_idx1);
                    if (bc.types_.contains(region.name()) and
                        bc.types_.at(region.name()) != BCType::zero_neumann)
                    {
                        if (bc.types_.at(region.name()) == BCType::dirichlet)
                        {
                            phi_facet = bc.values_.at(facet_idx);
                        }
                        else if (bc.types_.at(region.name()) == BCType::neumann)
                        {
                            /// @todo
                        }
                    }

                    // Add gradient contributions for the adjacent cell
                    for (int j = 0; j < dim; j++)
                    {
                        grad(cell_idx1, j) += phi_facet * area * normal(j);
                    }
                }
                // Internal facets
                else
                {
                    // Compute facet value from adjacent cell values
                    const real_t phi1 = phi(cell_idx1);
                    const real_t phi2 = phi(cell_idx2);
                    const real_t phi_facet = V->compute_facet_value(facet_idx, phi1, phi2);

                    // Add gradient contributions for the adjacent cells
                    for (int j = 0; j < dim; j++)
                    {
                        grad(cell_idx1, j) += phi_facet * area * normal(j);
                        grad(cell_idx2, j) -= phi_facet * area * normal(j);
                    }
                }
            }
        }

        // Gradient was incrementally constructed - needs assembly
        grad.assemble();

        // Scale by the (inverse) cell volume
        for (int cell_idx = 0; cell_idx < grad.n_owned(); cell_idx++)
        {
            const real_t vol_inv = 1.0 / V->cell_volume(cell_idx);
            for (int i = 0; i < grad.n_comp(); i++)
            {
                grad(cell_idx, i) *= vol_inv;
            }
        }

        // Ensure that ghost cell gradients are also updated
        grad.update_ghosts();
    }
    //=============================================================================
    void least_squares_gradient(const FVFunction &phi, FVFunction &grad)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const auto topo = mesh->topology();
        const int dim = mesh->pdim();

        // Zero the gradient
        // grad.set_all(0.0, true);
        grad.set_all(0.0);

        // Least-squares problem LHS matrix and RHS vector
        la::DenseMatrix A(dim, dim);
        la::DenseMatrix b(dim, 1);

        // Loop over all cells
        for (const auto &region : mesh->regions())
        {
            // Skip boundary regions
            if (region.dim() < dim)
            {
                continue;
            }

            for (const auto &[cell, cell_idx] : mesh->region_cells(region.name()))
            {
                // Reset least-squares LHS and RHS
                A.set_all(0.0);
                b.set_all(0.0);

                // Cell value and midpoint
                const real_t phi1 = phi(cell_idx);
                const auto r1 = V->cell_midpoint(cell_idx);

                for (const auto neighour : topo->adjacent_entities(cell_idx, dim, dim))
                {
                    // Neighbour value and midpoint
                    const real_t phi2 = phi(neighour);
                    const auto r2 = V->cell_midpoint(neighour);

                    // Value difference and intercell distance vector
                    const real_t dphi = phi2 - phi1;
                    const geo::Vec3 dr(r1, r2);

                    // Weighting factor
                    const real_t w = 1.0 / dr.mag();

                    // Add contributions for this neighbour
                    for (int i = 0; i < dim; i++)
                    {
                        for (int j = 0; j < dim; j++)
                        {
                            A(i, j) += w * dr(i) * dr(j);
                        }

                        b(i, 0) += w * dr(i) * dphi;
                    }
                }

                // Solve least squares problem and compute cell gradient
                const auto grad_cell = A.invert().first * b;
                for (int i = 0; i < dim; i++)
                {
                    grad(cell_idx, i) = grad_cell(i, 0);
                }
            }
        }
    }
    //=============================================================================
    void gradient(const FVFunction &phi, const FVBC &bc,
                  FVFunction &grad, GradientMethod method)
    {
        // Check that sizes match
        const real_t dim = phi.space()->mesh()->pdim();
        SFEM_CHECK_SIZES(phi.n_comp() * dim, grad.n_comp());

        switch (method)
        {
        case GradientMethod::green_gauss:
            green_gauss_gradient(phi, bc, grad);
            break;
        case GradientMethod::least_squares:
            least_squares_gradient(phi, grad);
            break;
        default:
            SFEM_ERROR(std::format("Invalid gradient computation method: {}\n",
                                   static_cast<int>(method)));
        }
    }
    //=============================================================================
    std::shared_ptr<FVFunction> gradient(const FVFunction &phi, const FVBC &bc,
                                         GradientMethod method)
    {
        // Quick access
        const auto V = phi.space();
        const auto mesh = V->mesh();
        const int dim = mesh->pdim();
        const int n_comp = phi.n_comp();

        // Set component names for the gradient
        std::vector<std::string> components(n_comp * dim);
        std::array<std::string, 3> suffixes = {"_x", "_y", "_z"};
        for (int i = 0; i < n_comp; i++)
        {
            const auto comp = phi.components()[i];
            for (int j = 0; j < dim; j++)
            {
                components[i * dim + j] = comp + suffixes[j];
            }
        }

        // Create, compute and return the gradient
        auto grad = std::make_shared<FVFunction>(V, components);
        gradient(phi, bc, *grad, method);
        return grad;
    }
}