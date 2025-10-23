#include "sfem.hpp"

using namespace sfem;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Read the mesh from file
    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);

    // Finite volume space and functions
    auto V = std::make_shared<fvm::FVSpace>(mesh);
    auto T = std::make_shared<fvm::FVField>(V, std::vector<std::string>{"T"});
    auto gradT = std::make_shared<fvm::FVField>(V, std::vector<std::string>{"Tx", "Ty"});

    // LHS matrix
    auto A = fvm::create_mat(*T);

    // RHS vector
    auto b = fvm::create_vec(*T);

    // Diffusion coefficient
    ConstantCoefficient coeff(1.0);

    std::vector<real_t> _vel = {0.0, 1.5};
    ConstantCoefficient vel(_vel);

    // Setup Dirichlet B.C.
    fvm::FVBC bc(T);
    bc.set_value("Left", "T", fvm::BCType::dirichlet, {1});
    bc.set_value("Right", "T", fvm::BCType::dirichlet, {0});

    // Solve the linear system using CG
    const real_t tol = 1e-5;
    const int n_iter_max = 2500;
    const bool verbose = true;
    la::CG solver(tol, n_iter_max, verbose);

    const int n_orthogonal_correctors = 1;
    for (int i = 0; i < n_orthogonal_correctors; i++)
    {
        log_msg(std::format("Non-Orthogonal Corrector - Iteration: {}\n", i), true);

        // Reset LHS and RHS
        A.set_all(0.0);
        b.set_all(0.0);

        // Assemble the LHS matrix and RHS vector
        fvm::laplacian(*T, *gradT, bc, coeff,
                       la::create_matset(A),
                       la::create_vecset(b));
        A.assemble();
        b.assemble();

        // Solve and update ghosted cell values
        solver.run(A, b, *T);
        T->update_ghosts();

        // Compute the gradient
        fvm::gradient(*T, bc, *gradT, fvm::GradientMethod::least_squares);

        // Save solution to file
        io::vtk::write(std::format("post/solution_00{}", i), *mesh, {T, gradT});
    }

    return 0;
}