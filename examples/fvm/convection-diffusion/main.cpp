// Solve the convection-diffusion equation on a 2D or 3D domain
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

    // Convective velocity
    ConstantCoefficient vel({0.0, 0.0});

    // Initial condition
    T->set_all(0.0);
    io::vtk::write("post/solution_0", *mesh, {T});

    // Boundary conditions
    fvm::FVBC bc(T);
    bc.set_value("Left", "T", fvm::BCType::dirichlet, {0});
    bc.set_value("Right", "T", fvm::BCType::dirichlet, {0});
    bc.set_value("Top", "T", fvm::BCType::dirichlet, {0});
    bc.set_value("Bottom", "T", fvm::BCType::dirichlet, {0});

    // LHS matrix
    auto A = fvm::create_mat(*T);

    // RHS vector
    auto b = fvm::create_vec(*T);

    // Linear solver
    la::SolverOptions options;
    options.tol = 1e-8;
    options.n_iter_max = 500;
    options.print_conv = true;
    la::GMRES solver(options);

    // Timestepping
    const real_t t_final = 50.0;
    const real_t t_start = 0.0;
    const real_t dt = 0.1;
    real_t time = t_start;
    int timestep = 0;
    while (time < t_final)
    {
        // Increment time
        time += dt;
        timestep++;
        log_msg(std::format("Time: {}, Timestep: {}\n", time, timestep), true);

        // Reset LHS and RHS
        A.set_all(0.0);
        b.set_all(0.0);

        // Add contributions
        fvm::laplacian(*T, *gradT, bc, nu,
                       la::create_matset(A),
                       la::create_vecset(b));
        fvm::convection(*T, bc, vel,
                        la::create_matset(A),
                        la::create_vecset(b));
        fvm::ddt(*T, ConstantCoefficient(1.0), dt,
                 la::create_matset(A),
                 la::create_vecset(b));

        // Assemble the LHS matrix and RHS vector
        A.assemble();
        b.assemble();

        // Solve and update ghosted cell values
        solver.run(A, b, *T);
        T->update_ghosts();

        // Save solution to file
        io::vtk::write(std::format("post/solution_{}", timestep), *mesh, {T});
    }

    return 0;
}