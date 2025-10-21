#include "fv_gradient.hpp"
#include "../../geo/utils.hpp"
#include "../../mesh/utils/geo_utils.hpp"
#include "../../mesh/utils/loop_utils.hpp"
#include "../../la/native/dense_matrix.hpp"

namespace sfem::fvm
{
    //=============================================================================
    void green_gauss_gradient(const FVField &phi,
                              const FVBC &bc,
                              FVField &grad)
    {
        // Zero the gradient
        grad.set_all(0.0);

        // Finite volume space
        const auto V = phi.space();

        // Construct gradient
        auto facet_work = [&](const mesh::Mesh &mesh,
                              const mesh::Region &region,
                              const mesh::Cell &,
                              int facet_idx)
        {
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
                        phi_facet = bc.values_.at(facet_idx)[0];
                    }
                    else if (bc.types_.at(region.name()) == BCType::neumann)
                    {
                        const real_t d = V->facet_cell_distances(facet_idx)[0];
                        phi_facet += bc.values_.at(facet_idx)[0] * d;
                    }
                    else if (bc.types_.at(region.name()) == BCType::robin)
                    {
                        /// @todo
                    }
                }

                // Add gradient contributions for the adjacent cell
                for (int j = 0; j < mesh.pdim(); j++)
                {
                    grad(cell_idx1, j) += phi_facet * area * normal(j);
                }
            }
            // Internal facets
            else
            {
                // Facet geometric interpolation factor
                const real_t g = V->facet_interp_factor(facet_idx);

                // Compute facet value from adjacent cell values
                const real_t phi1 = phi(cell_idx1);
                const real_t phi2 = phi(cell_idx2);
                const real_t phi_facet = g * phi1 + (1 - g) * phi2;

                // Add gradient contributions for the adjacent cells
                for (int j = 0; j < mesh.pdim(); j++)
                {
                    grad(cell_idx1, j) += phi_facet * area * normal(j);
                    grad(cell_idx2, j) -= phi_facet * area * normal(j);
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
            for (int i = 0; i < grad.n_comp(); i++)
            {
                grad(cell_idx, i) *= vol_inv;
            }
        };
        mesh::utils::for_all_cells(*V->mesh(), cell_work);

        // Ensure that ghost cell gradients are also updated
        grad.update_ghosts();
    }
    //=============================================================================
    void least_squares_gradient(const FVField &phi, FVField &grad)
    {
        // Zero the gradient
        grad.set_all(0.0);

        // Finite volume space, mesh and topology
        const auto V = phi.space();
        const auto mesh_ = V->mesh();
        const auto topo = mesh_->topology();

        // Least-squares problem LHS matrix and RHS vector
        la::DenseMatrix A(mesh_->pdim(), mesh_->pdim());
        la::DenseMatrix b(mesh_->pdim(), 1);

        // Construct gradient
        auto work = [&](const mesh::Mesh &mesh,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            // Reset least-squares LHS and RHS
            A.set_all(0.0);
            b.set_all(0.0);

            // Cell value and midpoint
            const real_t phi1 = phi(cell_idx);
            const auto r1 = V->cell_midpoint(cell_idx);

            for (const auto neighour : topo->adjacent_entities(cell_idx, mesh.pdim(), mesh.pdim()))
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
                for (int i = 0; i < mesh.pdim(); i++)
                {
                    for (int j = 0; j < mesh.pdim(); j++)
                    {
                        A(i, j) += w * dr(i) * dr(j);
                    }

                    b(i, 0) += w * dr(i) * dphi;
                }
            }

            // Solve least squares problem and compute cell gradient
            const auto grad_cell = A.invert().first * b;
            for (int i = 0; i < mesh.pdim(); i++)
            {
                grad(cell_idx, i) = grad_cell(i, 0);
            }
        };
        mesh::utils::for_all_cells(*mesh_, work);

        // Ensure that ghost cell gradients are also updated
        grad.update_ghosts();
    }
    //=============================================================================
    void gradient(const FVField &phi, const FVBC &bc,
                  FVField &grad, GradientMethod method)
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
    std::shared_ptr<FVField> gradient(const FVField &phi, const FVBC &bc,
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
        auto grad = std::make_shared<FVField>(V, components);
        gradient(phi, bc, *grad, method);
        return grad;
    }
}