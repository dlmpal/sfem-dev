#include <sfem/sfem.hpp>

using namespace sfem;
using namespace sfem::fem;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Mesh, order and FE space
    auto mesh = io::read_mesh(argv[1]);
    const int order = 2;
    const auto V = std::make_shared<fem::CGSpace>(mesh, order);

    // Temperature
    FEField T(V, {"T"});

    // Diffusivity
    ConstantField D("D", 1.0);

    // PETSc linear system
    auto Axb = fem::create_axb(T, la::SolverType::cg,
                               {}, la::Backend::petsc);

    // Create equation
    Equation eqn(T, Axb);
    eqn.add_kernel(Diffusion(T, D));

    // Set Dirichlet BC
    eqn.bc().set_value("Left", 10);
    eqn.bc().set_value("Right", 100);

    // Assemble and solve
    eqn.assemble();
    eqn.apply_dirichlet_bc();
    eqn.solve();

    // Write solution to file
    io::vtk::write("post/solution_000", {T});

    return 0;
}