#include "simple.hpp"

using namespace sfem;

void fill_bc(std::vector<fvm::FVBC> &Ubc, fvm::FVBC &Pbc)
{
    Ubc[0].set_value("Left", "Ux", fvm::BCType::dirichlet, {1.0});
    Ubc[0].set_value("Right", "Ux", fvm::BCType::zero_neumann, {});
    Ubc[0].set_value("Bottom", "Ux", fvm::BCType::dirichlet, {0.0});
    Ubc[0].set_value("Top", "Ux", fvm::BCType::dirichlet, {0.0});

    Ubc[1].set_value("Left", "Uy", fvm::BCType::dirichlet, {0.0});
    Ubc[1].set_value("Right", "Uy", fvm::BCType::zero_neumann, {});
    Ubc[1].set_value("Bottom", "Uy", fvm::BCType::dirichlet, {0.0});
    Ubc[1].set_value("Top", "Uy", fvm::BCType::dirichlet, {0.0});

    Pbc.set_value("Left", "P", fvm::BCType::zero_neumann, {});
    Pbc.set_value("Right", "P", fvm::BCType::dirichlet, {0.0});
    Pbc.set_value("Bottom", "P", fvm::BCType::zero_neumann, {});
    Pbc.set_value("Top", "P", fvm::BCType::zero_neumann, {});
}

void fill_ic(std::vector<std::shared_ptr<fvm::FVField>> &U, std::shared_ptr<fvm::FVField> P)
{
    for (auto u : U)
    {
        u->set_all(0.0);
    }
    P->set_all(0.0);
}

namespace sfem::la
{
    void jacobi_smoother(const SparseMatrix &A, const Vector &b, Vector &x)
    {
        Vector r(x.index_map(), x.block_size());
        spmv(A, x, r);  // r = Ax
        axpy(-1, b, r); // r -= b

        Vector d(x.index_map(), x.block_size());
        A.diagonal(d);

        // x += r / d
        for (int i = 0; i < r.index_map()->n_owned(); i++)
        {
            for (int j = 0; j < r.block_size(); j++)
            {
                x(i, j) -= r(i, j) / d(i, j);
            }
        }
    }
}

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Read the mesh from file
    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);

    // Finite volume space
    auto V = std::make_shared<fvm::FVSpace>(mesh);

    // Field names
    const std::vector<std::string> U_names = {"Ux", "Uy", "Uz"};
    const std::string P_name = "P";

    // Create fields and gradients
    std::vector<std::shared_ptr<fvm::FVField>> U(mesh->pdim());
    std::vector<std::shared_ptr<fvm::FVField>> U_old(mesh->pdim());
    std::vector<std::shared_ptr<fvm::FVField>> gradU(mesh->pdim());
    for (int i = 0; i < mesh->pdim(); i++)
    {
        U[i] = fvm::create_field(V, {U_names[i]});
        U_old[i] = fvm::create_field(V, {U_names[i]});
        gradU[i] = fvm::create_gradient(*U[i]);
    }
    auto P = fvm::create_field(V, {P_name});
    auto P_old = fvm::create_field(V, {P_name});
    auto gradP = fvm::create_gradient(*P);

    // Boundary conditions
    std::vector<fvm::FVBC> Ubc;
    for (int i = 0; i < mesh->pdim(); i++)
    {
        Ubc.emplace_back(U[i]);
    }
    fvm::FVBC Pbc(P);
    fill_bc(Ubc, Pbc);

    // Initial conditions
    fill_ic(U, P);

    // Aggregate the fields for plotting
    std::vector<std::shared_ptr<const Field>> fields_plot;
    for (int i = 0; i < mesh->pdim(); i++)
    {
        fields_plot.push_back(U[i]);
    }
    fields_plot.push_back(P);
    io::vtk::write("post/solution_0", *mesh, fields_plot);

    // Density
    const real_t rho = 1.0;
    const real_t rho_inv = 1.0 / rho;

    // Kinematic viscosity
    ConstantCoefficient nu(0.01);

    // Convective velocity
    FieldCoefficient vel(V->index_map(), {"Ux", "Uy"});

    // Pressure equation diffusivity coefficient
    FieldCoefficient D(V->index_map(), {"D"});

    // Momentum matrix diagonal and off-diagonal coefficients
    la::Vector a(V->index_map(), mesh->pdim());
    la::Vector H(V->index_map(), mesh->pdim());

    // LHS matrix and RHS vector
    // (Utilized for both momentum and pressure equations)
    auto A = fvm::create_mat(*P);
    auto b = fvm::create_vec(*P);

    // Momentum equation linear solver
    la::SolverOptions options_mom;
    options_mom.tol = 1e-6;
    options_mom.n_iter_max = 20;
    options_mom.print_conv = true;
    la::GMRES solver_mom(options_mom);

    // Pressure equation linear solver
    la::SolverOptions options_pres;
    options_pres.n_iter_max = 1000;
    options_pres.tol = 1e-8;
    options_pres.print_conv = true;
    la::CG solver_pres(options_pres);

    // Relaxation factors
    const real_t relax_factor_mom = 0.7;
    const real_t relax_factor_pres = 0.3;

    // SIMPLE configuration
    const int simple_max_iter = 60;
    const real_t simple_tol = 1e-8;
    int simple_iter = 0;

    // Plot interval
    const int plot_int = 10;

    // Timestepping
    const bool transient = false;
    const real_t dt_sgn = transient ? 1.0 : -1.0;
    const real_t dt = 0.01 * dt_sgn;
    const real_t t_final = 100.0;
    const real_t t_start = 0.0;
    real_t time = t_start;
    int timestep = 0;
    // Increment time
    // time += dt;
    // timestep++;
    // log_msg(std::format("Time: {}, Timestep: {}\n", time, timestep), true);

    while (simple_iter < simple_max_iter)
    {
        // Update fields
        for (int i = 0; i < mesh->pdim(); i++)
        {
            la::copy(*U[i], *U_old[i]);
        }
        la::copy(*P, *P_old);

        // Solve momentum equation
        for (int i = 0; i < mesh->pdim(); i++)
        {
            A.set_all(0.0);
            b.set_all(0.0);
            fvm::simple::assemble_momentum(*U[i], *gradU[i], Ubc[i], *gradP,
                                           nu, vel, rho, dt, i,
                                           la::create_matset(A),
                                           la::create_vecset(b));
            A.assemble();
            b.assemble();
            solver_mom.run(A, b, *U[i]);
            U[i]->update_ghosts();

            // Store diagonal
            A.diagonal(a, 0, i);
        }
        a.update_ghosts();

        /// @todo Make relaxation implicit
        // Relax momentum
        for (int i = 0; i < mesh->pdim(); i++)
        {
            fvm::simple::apply_relaxation(*U_old[i], *U[i], relax_factor_mom);
            U[i]->update_ghosts();
            fvm::gradient(*U[i], Ubc[i], *gradU[i], fvm::GradientMethod::least_squares);
        }

        // Assemble pressure diffusivity and RHS
        for (int i = 0; i < V->index_map()->n_local(); i++)
        {
            real_t a_sum = 0.0;
            for (int j = 0; j < mesh->pdim(); j++)
            {
                a_sum += a(i, j);
                H(i, j) = a(i, j) * (*U[j])(i) + rho_inv * (*gradP)(i, j) * V->cell_volume(j);
            }
            D(i, 0) = static_cast<real_t>(mesh->pdim()) / (a_sum);
        }

        // Solve pressure equation
        A.set_all(0.0);
        b.set_all(0.0);
        fvm::simple::assemble_pressure(*P, *gradP, Pbc, D, a, H,
                                       la::create_matset(A),
                                       la::create_vecset(b));
        A.assemble();
        b.assemble();
        for (int i = 0; i < 5; i++)
        {
            la::jacobi_smoother(A, b, *P);
        }
        solver_pres.run(A, b, *P);
        P->update_ghosts();

        // Relax pressure
        fvm::simple::apply_relaxation(*P_old, *P, relax_factor_pres);
        P->update_ghosts();

        // Compute pressure gradient
        fvm::gradient(*P, Pbc, *gradP, fvm::GradientMethod::least_squares);

        // Update convecting velocity (Rhie-Chow)
        for (int i = 0; i < V->index_map()->n_local(); i++)
        {
            for (int j = 0; j < mesh->pdim(); j++)
            {
                const real_t corr = 1 / a(i, j) * (H(i, j) - rho_inv * (*gradP)(i, j) * V->cell_volume(i));
                vel(i, j) = corr;
            }
        }

        // Save solution to file
        if (simple_iter % plot_int == 0)
        {
            io::vtk::write(std::format("post/solution_{}", simple_iter), *mesh, fields_plot);
        }

        // Check for convergence based on mass residual
        const real_t mass_res = fvm::simple::compute_mass_residual(gradU);
        log_msg(std::format("SIMPLE Iteration: {}, Mass residual: {}\n", simple_iter, mass_res), true);
        if (mass_res < simple_tol)
        {
            log_msg(std::format("SIMPLE converged at {} iterations. Mass residual {}\n", simple_iter, mass_res), true);
            break;
        }

        // Increment SIMPLE iteration
        simple_iter++;
    }

    return 0;
}