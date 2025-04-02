#include "sfem.hpp"
#include <iostream>

using namespace sfem;

namespace sfem::fem
{
    la::DenseMatrix plane_elasticity_kernel(const FEData &data,
                                            real_t E, real_t nu, real_t thick)
    {
        // Stress-strain matrix
        la::DenseMatrix D(3, 3);
        real_t coeff = thick * E / (1 - nu * nu);
        D(0, 0) = coeff * 1.0;
        D(0, 1) = coeff * nu;
        D(1, 0) = coeff * nu;
        D(1, 1) = coeff * 1.0;
        D(2, 2) = coeff * (1 - nu) * 0.5;

        // Strain-displacement matrix
        const auto &dNdX = data.dNdX;
        int n_cols = dNdX.n_rows() * 2;
        int n_rows = 3;
        la::DenseMatrix B(n_rows, n_cols);
        for (int i = 0; i < dNdX.n_rows(); i++)
        {
            // exx
            B(0, i * 2 + 0) = dNdX(i, 0);

            // eyy
            B(1, i * 2 + 1) = dNdX(i, 1);
            // exy
            B(2, i * 2 + 0) = dNdX(i, 1);
            B(2, i * 2 + 1) = dNdX(i, 0);
        }

        return B.T() * D * B;
    }
}

int main(int argc, char **argv)
{
    auto &app = initialize(argc, argv);

    auto mesh = io::read_mesh(argv[1]);
    int dim = mesh->topology()->dim();

    int order = 3;
    fem::CGSpace phi(mesh, order, {"Ux", "Uy"});

    real_t E = 5e9;
    real_t nu = 0.35;
    real_t thick = 0.5 * 1e-3;
    real_t pressure_value = 1000;

    // LHS matrix
    auto A = fem::petsc::create_mat(phi);
    for (int i = 0; i < mesh->topology()->n_entities(dim); i++)
    {
        // Only integrate locally owned cells
        if (mesh->topology()->entity_index_map(dim)->get_owner(i) != mpi::rank())
        {
            continue;
        }

        // Cell info
        auto cell = mesh->topology()->entity(i, dim);
        auto cell_dof = phi.index_map()->local_to_global(phi.cell_dof(i));
        auto cell_points = phi.cell_dof_points(i);

        // Integration
        auto element = phi.element(cell.type);
        la::DenseMatrix A_cell(cell_dof.size() * phi.n_comp(), cell_dof.size() * phi.n_comp());
        for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
        {
            auto data = element->transform(dim, element->integration_rule()->point(nqpt),
                                           cell_points);
            A_cell += fem::plane_elasticity_kernel(data, E, nu, thick) * data.detJ;
        }
        A.set_values(cell_dof, cell_dof, A_cell.data());
    }
    A.assemble();

    // RHS vector
    auto b = fem::petsc::create_vec(phi);
    for (auto region : {"Left"})
    {
        for (auto [facet, i] : mesh->region_facets(region))
        {
            // Only integrate locally owned facets
            if (mesh->topology()->entity_index_map(dim - 1)->get_owner(i) != mpi::rank())
            {
                continue;
            }

            // Facet info
            auto facet_points = phi.facet_dof_points(i);
            auto facet_normal = mesh::facet_normal(facet.type, facet_points).normalize();
            auto facet_dof = phi.index_map()->local_to_global(phi.facet_dof(i));

            // Integration
            auto element = phi.element(facet.type);
            la::DenseMatrix b_facet(facet_dof.size() * phi.n_comp(), 1);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                auto data = element->transform(dim,
                                               element->integration_rule()->point(nqpt),
                                               facet_points);

                for (int i = 0; i < data.N.n_rows(); i++)
                {
                    b_facet(i * phi.n_comp() + 0, 0) += -pressure_value * thick * data.N(i, 0) * data.detJ * facet_normal.x();
                    b_facet(i * phi.n_comp() + 1, 0) += -pressure_value * thick * data.N(i, 0) * data.detJ * facet_normal.y();
                }
            }
            b.set_values(facet_dof, b_facet.data());
        }
    }
    b.assemble();

    // Dirichlet B.C.
    fem::DirichletBC bc(phi);
    bc.set_value("Fixed", "Ux", 0);
    bc.set_value("Fixed", "Uy", 0);

    // Solution vector
    auto x = fem::petsc::create_vec(phi);

    // Apply Dirichlet BC and solve
    fem::petsc::solve(A, b, x, bc);

    // Save solution to VTK file
    io::vtk::write(std::format("post/solution_{}.vtk", mpi::rank()), phi, x.get_values());

    return 0;
}