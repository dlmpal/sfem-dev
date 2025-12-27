#include <sfem/sfem.hpp>

using namespace sfem;
using namespace sfem::fvm;

extern real_t conductivity(const std::array<real_t, 3> pt, real_t T);
extern real_t specific_heat_density(const std::array<real_t, 3> pt, real_t T);
extern void heat_source(const std::array<real_t, 3> pt, std::span<real_t> values, real_t time);

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Read the mesh from file
    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);

    // Finite volume space
    auto V = std::make_shared<FVSpace>(mesh);

    // Temperature field
    FVField T(V, {"T"});

    // Thermal conductivity
    FVField kappa(V, {"k"});

    // Density * heat capacity
    FVField rhocp(V, {"rhocp"});

    // Initial condition
    const real_t T_init = 25;
    T.values().set_all(T_init);

    // Boundary conditions
    const real_t h_inf = 30;
    const real_t T_inf = 25;
    const BCType bc_type = BCType::robin;
    const BCData bc_data{.a = h_inf, .b = 1.0, .c = T_inf};
    T.boundary_condition().set_region_bc("Left", bc_type, bc_data);
    T.boundary_condition().set_region_bc("Right", bc_type, bc_data);
    T.boundary_condition().set_region_bc("Front", bc_type, bc_data);
    T.boundary_condition().set_region_bc("Back", bc_type, bc_data);
    T.boundary_condition().set_region_bc("Top", bc_type, bc_data);
    T.boundary_condition().set_region_bc("Bottom", bc_type, bc_data);

    // Timestepping
    const real_t t_final = 500;
    const real_t t_start = 0.0;
    real_t dt = 0.1;
    real_t time = 0;
    int timestep = 0;
    bool after_weld_end = false;

    // Equation
    Equation eqn(T, create_axb(T, la::SolverType::cg));
    eqn.add_kernel(Laplacian(T, kappa))
        .add_kernel(ImplicitEuler(T, rhocp, dt))
        .add_kernel(Source(T, [&time](const FVField &T, int cell_idx, std::span<real_t> Q)
                           { heat_source(T.space()->cell_midpoint(cell_idx), Q, time); }));

    while (time <= t_final)
    {
        log_msg(std::format("Time: {}, Timestep: {}\n", time, timestep), true);

        // Update the coefficients
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            const auto x = V->cell_midpoint(cell_idx);
            kappa.cell_value(cell_idx) = conductivity(x, T.cell_value(cell_idx));
            rhocp.cell_value(cell_idx) = specific_heat_density(x, T.cell_value(cell_idx));
        };
        mesh::utils::for_all_cells(*mesh, work);

        // Solve
        eqn.assemble();
        eqn.solve();

        // Save solution to file
        io::vtk::write(std::format("post/solution_{}", timestep), {T});

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

real_t conductivity(const std::array<real_t, 3> pt, real_t T)
{
    real_t k = 0;
    if (T >= 20 and T < 800)
    {
        k = 54 - 3.33e-2 * T;
    }
    else if (T >= 800)
    {
        k = 27.3;
    }
    return k;
}

real_t specific_heat_density(const std::array<real_t, 3> pt, real_t T)
{
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
    return rho * cp;
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