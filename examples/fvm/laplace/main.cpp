#include "sfem.hpp"
#include <iostream>
#include <numeric>

using namespace sfem;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Read the mesh from file
    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);

    // Finite volume space and functions
    auto V = std::make_shared<fvm::FVSpace>(mesh);
    auto T = std::make_shared<fvm::FVFunction>(V, std::vector<std::string>{"T"});
    auto gradT = std::make_shared<fvm::FVFunction>(V, std::vector<std::string>{"Tx", "Ty"});

    // LHS matrix
    auto A = fvm::create_mat(*T);

    // RHS vector
    auto b = fvm::create_vec(*T);

    // Setup Dirichlet B.C.
    fvm::FVBC bc(T);
    bc.set_value("Left", "T", fvm::BCType::dirichlet, 10);
    bc.set_value("Right", "T", fvm::BCType::dirichlet, 100);

    // Assemble the LHS matrix and RHS vector
    ConstantCoefficient coeff(1.0);
    fvm::diffusion(*T, *gradT, bc, coeff,
                   la::create_matset(A),
                   la::create_vecset(b));
    A.assemble();
    b.assemble();
    b.update_ghosts();

    // Solve the resulting linear system using GMRES
    const real_t tol = 1e-5;
    const int n_iter_max = 600;
    const bool verbose = true;
    const int n_restart = 300;
    la::GMRES solver(tol, n_iter_max, verbose, n_restart);
    solver.run(A, b, *T);

    // Update ghosted cell values (required for plotting)
    T->update_ghosts();

    // Save solution to file
    io::vtk::write(std::format("post/solution_000"), *mesh, {T});

    return 0;
}