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
    auto phi = std::make_shared<fem::CGSpace>(mesh, order, std::vector<std::string>{"T"});

    // LHS matrix
    auto A = fem::petsc::create_mat(*phi);
    fem::assemble_matrix_cells(*phi, "Solid", fem::kernels::Diffusion3D(1.0), fem::petsc::create_matset(A));
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

    // Write solution to file
    io::vtk::write("post/solution_000", *phi, x.get_values());

    return 0;
}