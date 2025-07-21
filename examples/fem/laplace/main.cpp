#include "sfem.hpp"
#include <iostream>

using namespace sfem;

int main(int argc, char **argv)
{
    auto &app = initialize(argc, argv);
    // app.set_log_level(LogLevel::debug);

    // Mesh, order and FE space
    auto mesh = io::read_mesh(argv[1]);
    int order = 1;
    auto V = std::make_shared<fem::CGSpace>(mesh, order);

    // Temperature function
    auto T = std::make_shared<fem::FEFunction>(V, std::vector<std::string>{"T"});

    // Diffusion coefficient
    auto coeff = std::make_shared<ConstantCoefficient>(1.0);

    // LHS matrix
    auto A = fem::petsc::create_mat(*T);
    fem::assemble_matrix_cells(*T, "Solid", fem::kernels::Diffusion3D(coeff), la::petsc::create_matset(A));
    A.assemble();

    // RHS vector
    auto b = fem::petsc::create_vec(*T);
    b.assemble();

    // Dirichlet B.C.
    fem::DirichletBC bc(T);
    bc.set_value("Left", "T", 10);
    bc.set_value("Right", "T", 100);

    // Solution vector
    auto x = la::petsc::create_vec(*T);

    // Apply Dirichlet BC and solve
    fem::petsc::solve(A, b, x, bc);
    T->assemble();

    // Write solution to file
    io::vtk::write("post/solution_000", {T});

    return 0;
}