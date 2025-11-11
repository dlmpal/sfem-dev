#include "sfem.hpp"

using namespace sfem;

extern void conductivity(const std::array<real_t, 3> pt, std::span<real_t> values, real_t time);
extern void specific_heat_density(const std::array<real_t, 3> pt, std::span<real_t> values, real_t time);
extern void heat_source(const std::array<real_t, 3> pt, std::span<real_t> values, real_t time);

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

    // Thermal conductivity
    auto kappa = fvm::create_coeff(V, {"kappa"});

    // Density * heat capacity
    auto rhocp = fvm::create_coeff(V, {"rhocp"});

    // Initial condition
    const real_t T_init = 25;
    T->set_all(T_init);

    // Boundary conditions
    fvm::FVBC bc(T);
    const real_t h_inf = 30;
    const real_t T_inf = 25;
    bc.set_value("Left", "T", fvm::BCType::robin, {T_init, h_inf});
    bc.set_value("Right", "T", fvm::BCType::robin, {T_init, h_inf});
    bc.set_value("Front", "T", fvm::BCType::robin, {T_init, h_inf});
    bc.set_value("Back", "T", fvm::BCType::robin, {T_init, h_inf});
    bc.set_value("Top", "T", fvm::BCType::robin, {T_init, h_inf});
    bc.set_value("Bottom", "T", fvm::BCType::robin, {T_init, h_inf});

    // Source
    auto Q = create_field(V, {"Q"});

    // LHS matrix
    auto A = fvm::create_mat(*T);

    // RHS vector
    auto b = fvm::create_vec(*T);

    // Linear solver
    la::SolverOptions options;
    options.tol = 1e-8;
    options.n_iter_max = 500;
    options.print_conv = true;
    la::CG solver(options);

    // Timestepping
    const real_t t_final = 500;
    const real_t t_start = 0.0;
    real_t dt = 0.1;
    real_t time = 0;
    int timestep = 0;
    bool after_weld_end = false;
    while (time <= t_final)
    {
        log_msg(std::format("Time: {}, Timestep: {}\n", time, timestep), true);

        // Setup coefficients
        la::copy(*T, *kappa.field());
        la::copy(*T, *rhocp.field());
        fvm::eval_field(*kappa.fv_field(), conductivity, time);
        fvm::eval_field(*rhocp.fv_field(), specific_heat_density, time);

        // Reset LHS and RHS
        A.set_all(0.0);
        b.set_all(0.0);

        // Add diffusion contribution
        fvm::laplacian(*T, *gradT, bc, kappa,
                       la::create_matset(A),
                       la::create_vecset(b));
        fvm::ddt(*T, rhocp, dt,
                 la::create_matset(A),
                 la::create_vecset(b));

        // Add source term contribution
        fvm::eval_field(*Q, heat_source, time);
        fvm::add_source_term(*Q, la::create_vecset(b));

        // Assemble the LHS matrix and RHS vector
        A.assemble();
        b.assemble();

        // Solve and update ghost cell values
        solver.run(A, b, *T);
        T->update_ghosts();

        // Compute the gradient
        fvm::gradient(*T, bc, *gradT);

        // Save solution to file
        io::vtk::write(std::format("post/solution_{}", timestep), *mesh, {T, gradT});

        time += dt;
        timestep++;

        if (time > 25 and after_weld_end == false)
        {
            after_weld_end = true;
            dt *= 100;
        }
    }

    return 0;
}

void conductivity(const std::array<real_t, 3> pt, std::span<real_t> values, real_t time)
{
    const real_t T = values[0];
    real_t k = 0;
    if (T >= 20 and T < 800)
    {
        k = 54 - 3.33e-2 * T;
    }
    else if (T >= 800)
    {
        k = 27.3;
    }
    values[0] = k;
}

void specific_heat_density(const std::array<real_t, 3> pt, std::span<real_t> values, real_t time)
{
    const real_t T = values[0];
    const real_t rho = 8050;
    real_t cp = 0;
    if (T >= 20 and T < 600)
    {
        cp = 425 + (7.73e-1 - 1.69e-3 * T + 2.22e-6 * T * T) * T;
    }
    else if (T >= 600 and T < 735)
    {
        cp = 666 + 13002 / (738 - T);
    }
    else if (T >= 735 and T < 900)
    {
        cp = 545 + 17820 / (T - 731);
    }
    else if (T >= 900)
    {
        cp = 650;
    }
    values[0] = rho * cp;
}

void heat_source(const std::array<real_t, 3> pt, std::span<real_t> values, real_t time)
{
    // Constant parameters
    const real_t Q = 990;
    const real_t vel = 0.0033;
    const real_t center_y = 0.025;
    const real_t center_z = 0.0;
    const real_t sigma = 2e-3;
    const real_t t_weld = 0.04 / vel;

    // Compute distannce from electrode
    const real_t r = geo::compute_distance(pt, {vel * time, center_y, center_z});

    // Compute the source value
    const real_t coeff = (Q / vel) / (2 * sigma * sigma * M_PI);

    if (time < t_weld)
    {
        values[0] = coeff * std::exp(-r * r / (2 * sigma * sigma));
    }
    else
    {
        values[0] = 0;
    }
}