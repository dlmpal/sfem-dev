#include "fv_gradient.hpp"
#include <iostream>
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
        const int n_comp = phi.n_comp();

        // Zero the gradient
        grad.set_all(0.0);

        FVFunction grad_est(V, grad.components());
        grad_est.set_all(0.0);

        // Integrate interior facets
        for (const auto &region : mesh->regions())
        {
            for (const auto &[facet, facet_idx] : mesh->region_facets(region.name()))
            {
                // Integrate locally owned facets only
                if (mesh->topology()->entity_index_map(mesh->pdim() - 1)->is_ghost(facet_idx))
                {
                    continue;
                }

                // Facet's adjacent cells
                const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
                const auto &[cell_idx1, cell_idx2] = adjacent_cells;

                // Facet area and normal vector
                const real_t area = V->facet_area(facet_idx);
                const auto normal = V->facet_normal(facet_idx);

                // Boundary facets
                if (cell_idx1 == cell_idx2)
                {
                    for (int i = 0; i < n_comp; i++)
                    {
                        // Boundary facet value
                        // Corresponds to zero Neumann BC by default
                        real_t value = phi(cell_idx1, i);
                        if (bc.types_.contains(region.name()) and
                            bc.types_.at(region.name()) != BCType::zero_neumann)
                        {
                            if (bc.types_.at(region.name()) == BCType::dirichlet)
                            {
                                value = bc.values_.at(facet_idx * n_comp + i);
                            }
                            else if (bc.types_.at(region.name()) == BCType::neumann)
                            {
                                /// @todo
                            }
                        }

                        for (int j = 0; j < dim; j++)
                        {
                            grad_est(cell_idx1, j) += value * area * normal(j);
                        }
                    }
                }
                // Internal facets
                else
                {
                    for (int i = 0; i < n_comp; i++)
                    {
                        // Adjacent cell values and facet value
                        const real_t phi1 = phi(adjacent_cells[0], i);
                        const real_t phi2 = phi(adjacent_cells[1], i);
                        const real_t phi_facet = 0.5 * (phi1 + phi2);
                        // const real_t phi_facet = V->compute_facet_value(facet_idx,
                        // phi1, phi2);

                        for (int j = 0; j < dim; j++)
                        {
                            grad_est(cell_idx1, j) += phi_facet * area * normal(j);
                            grad_est(cell_idx2, j) -= phi_facet * area * normal(j);
                        }
                    }
                }
            }
        }

        // Scale by (inverse) cell volume
        for (int cell_idx = 0; cell_idx < grad.n_local(); cell_idx++)
        {
            const real_t vol_inv = 1 / V->cell_volume(cell_idx);
            for (int i = 0; i < grad.n_comp(); i++)
            {
                grad_est(cell_idx, i) *= vol_inv;
            }
        }

        // Integrate interior facets
        for (const auto &region : mesh->regions())
        {
            for (const auto &[facet, facet_idx] : mesh->region_facets(region.name()))
            {
                // Integrate locally owned facets only
                if (mesh->topology()->entity_index_map(mesh->pdim() - 1)->is_ghost(facet_idx))
                {
                    continue;
                }

                // Facet's adjacent cells
                const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
                const auto &[cell_idx1, cell_idx2] = adjacent_cells;

                // Facet area and normal vector
                const real_t area = V->facet_area(facet_idx);
                const auto normal = V->facet_normal(facet_idx);

                // Boundary facets
                if (cell_idx1 == cell_idx2)
                {
                    for (int i = 0; i < n_comp; i++)
                    {
                        // Boundary facet value
                        // Corresponds to zero Neumann BC by default
                        real_t value = phi(cell_idx1, i);
                        if (bc.types_.contains(region.name()) and
                            bc.types_.at(region.name()) != BCType::zero_neumann)
                        {
                            if (bc.types_.at(region.name()) == BCType::dirichlet)
                            {
                                value = bc.values_.at(facet_idx * n_comp + i);
                            }
                            else if (bc.types_.at(region.name()) == BCType::neumann)
                            {
                                /// @todo
                            }
                        }

                        for (int j = 0; j < dim; j++)
                        {
                            grad(cell_idx1, j) += value * area * normal(j);
                        }
                    }
                }
                // Internal facets
                else
                {
                    for (int i = 0; i < n_comp; i++)
                    {
                        // Adjacent cell values and facet value
                        const real_t phi1 = phi(adjacent_cells[0], i);
                        const real_t phi2 = phi(adjacent_cells[1], i);
                        const geo::Vec3 grad1(grad_est(cell_idx1, 0), grad_est(cell_idx1, 1), 0);
                        const geo::Vec3 grad2(grad_est(cell_idx2, 0), grad_est(cell_idx2, 1), 0);
                        geo::Vec3 rf({}, V->facet_midpoint(facet_idx));
                        geo::Vec3 r1({}, V->cell_midpoint(cell_idx1));
                        geo::Vec3 r2({}, V->cell_midpoint(cell_idx2));

                        const real_t cor = geo::inner(grad1 + grad2, rf - 0.5 * (r1 + r2));
                        const real_t phi_facet = 0.5 * (phi1 + phi2) + 0.5 * cor;

                        for (int j = 0; j < dim; j++)
                        {
                            grad(cell_idx1, j) += phi_facet * area * normal(j);
                            grad(cell_idx2, j) -= phi_facet * area * normal(j);
                        }
                    }
                }
            }
        }

        // Scale by (inverse) cell volume
        for (int cell_idx = 0; cell_idx < grad.n_local(); cell_idx++)
        {
            const real_t vol_inv = 1 / V->cell_volume(cell_idx);
            for (int i = 0; i < grad.n_comp(); i++)
            {
                grad(cell_idx, i) *= vol_inv;
            }
        }
    }
    //=============================================================================
    void gradient(const FVFunction &phi,
                  const FVBC &bc,
                  FVFunction &grad,
                  GradientMethod method)
    {
        // Check that sizes match
        const real_t dim = phi.space()->mesh()->pdim();
        SFEM_CHECK_SIZES(phi.n_comp() * dim, grad.n_comp());

        switch (method)
        {
        case GradientMethod::green_gauss:
            green_gauss_gradient(phi, bc, grad);
            break;
        default:
            SFEM_ERROR(std::format("Invalid gradient computation method: {}\n",
                                   static_cast<int>(method)));
        }
    }
    //=============================================================================
    std::shared_ptr<FVFunction> gradient(const FVFunction &phi,
                                         const FVBC &bc,
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