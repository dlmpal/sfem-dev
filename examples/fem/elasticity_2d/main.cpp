#include "sfem.hpp"
#include <iostream>

using namespace sfem;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Mesh, order and FE space
    auto mesh = io::read_mesh(argv[1]);
    const int order = 2;
    auto V = std::make_shared<fem::CGSpace>(mesh, order);

    // Displacement function
    auto U = std::make_shared<fem::FEFunction>(V, std::vector<std::string>{"Ux", "Uy"});

    // Elasticity coefficients and pressure
    auto thick = std::make_shared<ConstantCoefficient>(0.5 * 1e-3);
    auto E = std::make_shared<ConstantCoefficient>(5e9);
    auto nu = std::make_shared<ConstantCoefficient>(0.35);
    auto pressure = std::make_shared<ConstantCoefficient>(1000);

    // LHS matrix
    auto A = fem::petsc::create_mat(*U);
    fem::assemble_matrix_cells(*U, "Solid", fem::kernels::LinearElasticity2D(E, nu, thick), la::petsc::create_matset(A));
    A.assemble();

    // RHS vector
    auto b = fem::petsc::create_vec(*U);
    fem::assemble_vec_facets(*U, "Left", fem::kernels::PressureLoad2D(thick, pressure), la::petsc::create_vecset(b));
    b.assemble();

    // Dirichlet B.C.
    fem::DirichletBC bc(U);
    bc.set_value("Fixed", "Ux", 0);
    bc.set_value("Fixed", "Uy", 0);

    // Solution vector
    auto x = la::petsc::create_vec(*U);

    // Apply Dirichlet BC and solve
    fem::petsc::solve(A, b, x, bc);
    U->update_ghosts();

    // Save solution to VTK file
    io::vtk::write("post/solution_000", {U});

    return 0;
}