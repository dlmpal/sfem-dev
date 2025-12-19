#include <sfem/sfem.hpp>

using namespace sfem;
using namespace sfem::fem;
using namespace sfem::fem::solid_mechanics;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Mesh, order and FE space
    auto mesh = io::read_mesh(argv[1]);
    const int order = 2;
    auto V = std::make_shared<fem::CGSpace>(mesh, order);

    // Displacement field
    FEField U(V, {"Ux", "Uy"});

    // Elasticity coefficients and pressure
    ConstantField E("thick", 5e9);
    ConstantField nu("thick", 0.35);
    ConstantField rho("rho", 0.0);
    ConstantField P("P", 1e3);

    // Constitutive law
    LinearElasticPlaneStress constitutive(E, nu, rho);

    // Strain
    SmallStrain strain(U);

    // Stress
    Stress stress(U, strain, constitutive);

    // PETSc linear system
    auto Axb = fem::create_axb(U, la::SolverType::cg,
                               {.n_iter_max = 1000},
                               la::Backend::petsc);

    // Linear elasticity kernel
    LinearElasticity elasticity(U, strain, constitutive);
    PressureLoad pressure_load(U, P, mesh->get_region_by_name("Left"));

    // Create equation
    Equation eqn(U, Axb);
    eqn.add_kernel(elasticity);
    eqn.add_kernel(pressure_load);

    // Set Dirichlet BC
    eqn.bc().set_value("Fixed", 0, 0);
    eqn.bc().set_value("Fixed", 0, 1);

    // Assemble and solve
    eqn.assemble();
    eqn.apply_dirichlet_bc();
    eqn.solve();

    // Compute cell strains
    FEField Eps(std::make_shared<fem::CGSpace>(mesh, 0), {"exx", "eyy", "exy"});
    cell_qpoint_average(*V, strain, Eps);

    // Compute cell stresses
    FEField S(std::make_shared<fem::CGSpace>(mesh, 0), {"sxx", "syy", "sxy"});
    cell_qpoint_average(*V, stress, S);

    // Save solution to VTK file
    io::vtk::write("post/solution_000", {U}, {S, Eps});

    return 0;
}
