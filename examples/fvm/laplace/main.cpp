// Solve Laplace's equation on a 2D or 3D domain
// Create the mesh with the cart-mesh application by:
// ${SFEM_DEV_INSTALL_DIR}/bin/cart-mesh -d=2 -Nx=100 -Ny=100 -x-low=0 -x-high=1 -y-low=0 -y-high=1

#include "sfem.hpp"

using namespace sfem;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Read the mesh from file
    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);

    // Finite volume space
    auto V = std::make_shared<fvm::FVSpace>(mesh);

    // Temperature field and gradient
    auto T = fvm::create_field(V, {"T"});
    auto gradT = fvm::create_gradient(*T);

    // Diffusivity
    ConstantCoefficient nu(1.0);

    // LHS matrix
    auto A = fvm::create_mat(*T);

    // RHS vector
    auto b = fvm::create_vec(*T);

    // Boundary conditions
    fvm::FVBC bc(T);
    bc.set_value("Left", "T", fvm::BCType::dirichlet, {1});
    bc.set_value("Right", "T", fvm::BCType::dirichlet, {0});

    // Linear solver
    la::SolverOptions options;
    options.print_iter = true;
    la::CG solver(options);

    const int n_orthogonal_correctors = 2;
    for (int i = 0; i < n_orthogonal_correctors; i++)
    {
        log_msg(std::format("Non-Orthogonal Corrector - Iteration: {}\n", i), true);

        // Reset LHS and RHS
        A.set_all(0.0);
        b.set_all(0.0);

        // Assemble the LHS matrix and RHS vector
        fvm::laplacian(*T, *gradT, bc, nu,
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