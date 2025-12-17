#include "sfem.hpp"

using namespace sfem;
using namespace sfem::fem;
using namespace sfem::fem::kernels;
using namespace sfem::fem::kernels::elasticity;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Mesh, order and FE space
    auto mesh = io::read_mesh(argv[1]);
    const int order = 2;
    auto V = std::make_shared<fem::CGSpace>(mesh, order);

    // Displacement field
    FEField U(V, {"Ux", "Uy"});

    // Dirichlet B.C.
    fem::DirichletBC bc(V, U.n_comp());
    bc.set_value("Fixed", 0, 0);
    bc.set_value("Fixed", 0, 1);
    bc.set_value("Upper", 0.01, 1);

    // Elasticity coefficients and pressure
    ConstantField E("thick", 5e9);
    ConstantField nu("thick", 0.35);
    ConstantField rho("rho", 0.0);
    ConstantField thick("thick", 0.5 * 1e-3);
    // ConstantField P("P", 1e3);

    // Constitutive law
    LinearElasticPlaneStress constitutive(E, nu, rho, thick);

    auto Axb = fem::create_axb(U, la::SolverType::cg,
                               {.rtol = 1e-10, .n_iter_max = 1000},
                               la::Backend::petsc);

    LinearElasticity elasticity(U, constitutive);
    elasticity(Axb->lhs(), Axb->rhs());
    Axb->assemble();

    const auto [dof_idxs, dof_values] = bc.get_dofs_values();
    Axb->eliminate_dofs(dof_idxs, dof_values);

    Axb->solve(U.dof_values());
    U.dof_values().update_ghosts();

    // Save solution to VTK file
    io::vtk::write("post/solution_000", {U});

    return 0;
}