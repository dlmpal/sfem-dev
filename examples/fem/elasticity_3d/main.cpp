#include "sfem.hpp"
#include "argparse.hpp"
#include <iostream>

using namespace sfem;
using namespace sfem::fem;
using namespace sfem::fem::kernels;
using namespace sfem::fem::kernels::elasticity;
using namespace argparse;

extern void solve_static(FEField &U, Strain &strain, LinearElasticIsotropic &constitutive);
extern void solve_modal(const FEField &U, Strain &strain, LinearElasticIsotropic &constitutive, int n_eigs);

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Parse command line arguments
    ArgParser parser;
    parser.add_argument(Argument("mesh-dir", true));
    parser.add_argument(Argument("order", true));
    parser.add_argument(Argument("problem-type", true));
    parser.add_argument(Argument("n-eigs").value(6));
    parser.parse_args(argc, argv);

    // Mesh, order and FE space
    const auto mesh = io::read_mesh(parser.get_argument("mesh-dir")->value<std::string>());
    const int order = parser.get_argument("order")->value<int>();
    const auto V = std::make_shared<fem::CGSpace>(mesh, order);

    // Displacement field
    FEField U(V, {"Ux", "Uy", "Uz"});

    // Elasticity coefficients
    ConstantField E("E", 1e5);
    ConstantField nu("nu", 0.2);
    ConstantField rho("rho", 1e-3);

    // Constitutive law
    LinearElastic3D constitutive(E, nu, rho);

    // Strain
    Strain strain(U);

    // Get problem type and solve
    auto problem_type = parser.get_argument("problem-type")->value<std::string>();
    if (problem_type == "modal")
    {
        const int n_eigs = parser.get_argument("n-eigs")->value<int>();
        solve_modal(U, strain, constitutive, 6);
    }
    else if (problem_type == "static")
    {
        solve_static(U, strain, constitutive);
    }
    else
    {
        SFEM_ERROR(std::format("Invalid problem type: {}\n", problem_type));
    }

    return 0;
}

void set_bc(DirichletBC &bc)
{
    bc.set_value("Left", 0, 0);
    bc.set_value("Left", 0, 1);
    bc.set_value("Left", 0, 2);
    bc.set_value("Right", 0, 0);
    bc.set_value("Right", 0, 1);
    bc.set_value("Right", 0, 2);
}

void solve_static(FEField &U, Strain &strain, LinearElasticIsotropic &constitutive)
{
    // Quick access
    const auto V = U.space();

    LinearElasticity elasticity(U, strain, constitutive, {0, -9.81 * 5, 0});

    auto Axb = fem::create_axb(U, la::SolverType::cg,
                               {.rtol = 1e-10, .n_iter_max = 1500},
                               la::Backend::petsc);

    Equation eqn(U, Axb);
    eqn.add_kernel(elasticity);
    set_bc(eqn.bc());

    eqn.assemble();
    eqn.apply_dirichlet_bc();
    eqn.solve();

    Stress stress(U, strain, constitutive);

    FEField E(std::make_shared<CGSpace>(V->mesh(), 0), {"exx", "eyy", "ezz", "exy", "eyz", "exz"});
    cell_qpoint_average(*V, strain, E);
    FEField S(std::make_shared<CGSpace>(V->mesh(), 0), {"sxx", "syy", "szz", "sxy", "syz", "sxz"});
    cell_qpoint_average(*V, stress, S);

    // Save solution to VTK file
    io::vtk::write("post/solution_000", {U}, {E, S});
}

void solve_modal(const FEField &U, Strain &strain, LinearElasticIsotropic &constitutive, int n_eigs)
{
    // Mass matrix
    auto M = fem::petsc::create_mat(U);
    MassND mass(U, constitutive.rho());
    mass(la::petsc::create_matset(M), {});
    M.assemble();

    // Stiffness matrix
    auto K = fem::petsc::create_mat(U);
    LinearElasticity elasticity(U, strain, constitutive);
    elasticity(la::petsc::create_matset(K), {});
    K.assemble();

    // Apply Dirichlet BC
    DirichletBC bc(U.space(), U.n_comp());
    set_bc(bc);
    la::petsc::eliminate_rows_cols(bc.get_dofs_values().first, K);

    // Solve the eigenproblem with SLEPc
    la::slepc::SlepcEPS eps;
    eps.set_from_options();
    eps.set_operators(K, M);
    eps.solve(n_eigs);

    // Print the computed eigenfrequencies
    if (mpi::rank() == 0)
    {
        for (int i = 0; i < std::min(n_eigs, eps.n_converged()); i++)
        {
            auto [eigval_real, eigval_imag] = eps.eigenvalue(i);
            std::cout << std::format("Mode: {}, Frequency: {} [Hz]\n", i, std::sqrt(eigval_real) / 2 / M_PI);
        }
    }
}