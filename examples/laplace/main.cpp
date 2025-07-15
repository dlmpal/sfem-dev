#include "sfem.hpp"
#include <iostream>

using namespace sfem;

int main(int argc, char **argv)
{
    auto &app = initialize(argc, argv);
    // app.set_log_level(LogLevel::debug);

    auto mesh = io::read_mesh(argv[1]);
    int dim = mesh->topology()->dim();

    int order = 2;
    auto V = std::make_shared<fem::CGSpace>(mesh, order);
    auto T = std::make_shared<fem::FEFunction>(V, std::vector<std::string>{"T"});

    // Diffusion coefficient
    auto coeff = std::make_shared<fem::ConstantCoefficient>(1.0);

    // LHS matrix
    auto A = fem::petsc::create_mat(*V, 1);
    fem::assemble_matrix_cells(*V, "Solid", 1, fem::kernels::Diffusion3D(coeff), fem::petsc::create_matset(A));
    A.assemble();

    // RHS vector
    auto b = fem::petsc::create_vec(*V, 1);
    b.assemble();

    // Dirichlet B.C.
    fem::DirichletBC bc(T);
    bc.set_value("Fixed", "T", 10);
    bc.set_value("Free", "T", 100);

    // Solution vector
    auto x = la::petsc::create_vec(*T);

    // Apply Dirichlet BC and solve
    fem::petsc::solve(A, b, x, bc);
    T->assemble();

    // Write solution to file
    io::vtk::write("post/solution_000", {T});

    return 0;
}