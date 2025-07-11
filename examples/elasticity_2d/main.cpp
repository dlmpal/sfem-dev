#include "sfem.hpp"
#include <iostream>

using namespace sfem;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    auto mesh = io::read_mesh(argv[1]);
    int dim = mesh->topology()->dim();

    const int order = 1;
    auto phi = std::make_shared<fem::CGSpace>(mesh, order, std::vector<std::string>{"Ux", "Uy"});

    // Elasticity coefficients and pressure
    auto thick = std::make_shared<ConstantFunction>(0.5 * 1e-3);
    auto E = std::make_shared<ConstantFunction>(5e9);
    auto nu = std::make_shared<ConstantFunction>(0.35);
    auto pressure = std::make_shared<ConstantFunction>(1000);

    // LHS matrix
    auto A = fem::petsc::create_mat(*phi);
    fem::assemble_matrix_cells(*phi, "Solid",
                               fem::kernels::LinearElasticity2D(E, nu, thick), fem::petsc::create_matset(A));
    A.assemble();

    // RHS vector
    auto b = fem::petsc::create_vec(*phi);
    fem::assemble_vec_facets(*phi, "Left",
                             fem::kernels::PressureLoad2D(thick, pressure), fem::petsc::create_vecset(b));
    b.assemble();

    // Dirichlet B.C.
    fem::DirichletBC bc(phi);
    bc.set_value("Fixed", "Ux", 0);
    bc.set_value("Fixed", "Uy", 0);

    // Solution vector
    auto x = fem::petsc::create_vec(*phi);

    // Apply Dirichlet BC and solve
    fem::petsc::solve(A, b, x, bc);

    // Save solution to VTK file
    io::vtk::write("post/solution_000", *phi, x.get_values());

    return 0;
}