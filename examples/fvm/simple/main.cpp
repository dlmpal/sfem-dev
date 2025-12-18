// Solve the Navier-Stokes equations for 2D laminar channel flow
// Create the mesh with the cart-mesh application by:
// ${SFEM_DEV_INSTALL_DIR}/bin/cart-mesh -d=2 -Nx=200 -Ny=20 -x-low=0 -x-high=1 -y-low=0 -y-high=0.1

#include "simple.hpp"

extern void set_bc(std::vector<FVField> U, FVField &P);

int main(int argc, char *argv[])
{
    initialize(argc, argv, "SIMPLESolver");

    // Read mesh
    auto mesh = io::read_mesh("mesh", mesh::PartitionCriterion::shared_facet);

    // Finite volume space
    auto V = std::make_shared<FVSpace>(mesh);

    // Create fields and set BC
    auto [U, P] = create_fields(V);
    set_bc(U, P);

    // Density and dynamic viscosity
    const real_t rho = 1;
    const real_t mu = 1e3;

    // SIMPLE solver options
    algo::SIMPLEOptions options;
    options.pressure_solver_options.rtol = 1e-2;

    // Create the solver
    algo::SIMPLESolver solver(U, P, rho, mu, options);

    for (int iter = 0; iter < options.max_iter_simple; iter++)
    {
        log_msg(std::format("SIMPLE Iteration: {}\n", iter), true);

        solver.step(0.0, 0.0);

        if (iter % options.plot_int == 0)
        {
            auto fields_to_plot = U;
            fields_to_plot.push_back(P);
            io::vtk::write(std::format("post/solution_{}", iter), fields_to_plot);
        }
    }

    return 0;
}

void set_bc(std::vector<FVField> U, FVField &P)
{
    auto &ubc = U[0].boundary_condition();
    ubc.set_region_bc("Left", BCType::dirichlet, 0.01);
    ubc.set_region_bc("Right", BCType::zero_neumann, 0.0);
    ubc.set_region_bc("Bottom", BCType::dirichlet, 0.0);
    ubc.set_region_bc("Top", BCType::dirichlet, 0.0);

    auto &vbc = U[1].boundary_condition();
    vbc.set_region_bc("Left", BCType::dirichlet, 0.0);
    vbc.set_region_bc("Right", BCType::zero_neumann, 0.0);
    vbc.set_region_bc("Bottom", BCType::dirichlet, 0.0);
    vbc.set_region_bc("Top", BCType::dirichlet, 0.0);

    auto &pbc = P.boundary_condition();
    pbc.set_region_bc("Left", BCType::zero_neumann, 0.0);
    pbc.set_region_bc("Right", BCType::dirichlet, 0.0);
    pbc.set_region_bc("Bottom", BCType::zero_neumann, 0.0);
    pbc.set_region_bc("Top", BCType::zero_neumann, 0.0);
}