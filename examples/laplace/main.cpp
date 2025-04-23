#include "sfem.hpp"
#include <iostream>

using namespace sfem;

int main(int argc, char **argv)
{
    auto &app = initialize(argc, argv);
    // app.set_log_level(LogLevel::debug);

    auto mesh = io::read_mesh(argv[1]);
    int dim = mesh->topology()->dim();

    int order = 1;
    int n_vars = 1;
    auto phi = std::make_shared<fem::CGSpace>(mesh, n_vars, std::vector<std::string>{"T"});

    // LHS matrix
    auto A = fem::petsc::create_mat(*phi);
    for (int i = 0; i < mesh->topology()->n_entities(dim); i++)
    {
        // Only integrate locally owned cells
        if (mesh->topology()->entity_index_map(dim)->get_owner(i) != mpi::rank())
        {
            continue;
        }

        // Cell info
        auto cell = mesh->topology()->entity(i, dim);
        auto cell_dof = phi->index_map()->local_to_global(phi->cell_dof(i));
        auto cell_points = phi->cell_dof_points(i);

        // Integration
        auto element = phi->element(cell.type);
        la::DenseMatrix A_cell(cell_dof.size(), cell_dof.size());
        for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
        {
            auto data = element->transform(dim,
                                           element->integration_rule()->point(nqpt),
                                           cell_points);
            for (int i = 0; i < data.N.n_rows(); i++)
            {
                for (int j = 0; j < data.N.n_rows(); j++)
                {
                    for (int k = 0; k < dim; k++)
                    {
                        A_cell(i, j) += data.dNdX(i, k) * data.dNdX(j, k) * data.detJ;
                    }
                }
            }
        }
        A.set_values(cell_dof, cell_dof, A_cell.data());
    }
    A.assemble();

    // RHS vector
    auto b = fem::petsc::create_vec(*phi);
    b.assemble();

    // Dirichlet B.C.
    fem::DirichletBC bc(phi);
    bc.set_value("Left", "T", 10);
    bc.set_value("Right", "T", 100);

    // Solution vector
    auto x = fem::petsc::create_vec(*phi);

    // Apply Dirichlet BC and solve
    fem::petsc::solve(A, b, x, bc);

    io::vtk::write("post/solution_000", *phi, x.get_values());

    return 0;
}