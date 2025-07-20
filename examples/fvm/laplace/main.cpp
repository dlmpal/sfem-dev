#include "sfem.hpp"
#include <iostream>

using namespace sfem;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);
    auto topology = mesh->topology();

    auto V = std::make_shared<fvm::FVSpace>(mesh);
    auto T = std::make_shared<fvm::FVFunction>(V, std::vector<std::string>{"T"});

    // LHS matrix
    auto A = fvm::petsc::create_mat(*T);

    // RHS vector
    auto b = fvm::petsc::create_vec(*T);

    ConstantCoefficient coeff(1.0);

    fvm::isotropic_diffusion(*T, coeff, "Solid", la::petsc::create_matset(A), la::petsc::create_vecset(b));
    fvm::dirichlet_bc(*T, coeff, "Left", 10, la::petsc::create_matset(A), la::petsc::create_vecset(b));
    fvm::dirichlet_bc(*T, coeff, "Right", 100, la::petsc::create_matset(A), la::petsc::create_vecset(b));

    A.assemble();
    b.assemble();

    auto x = la::petsc::create_vec(*T);
    x.assemble();

    la::petsc::solve(A, b, x);
    T->assemble();

    io::vtk::write("post/solution_000", *mesh, {T});

    return 0;
}