#include "sfem.hpp"

using namespace sfem;
using namespace sfem::fem;

int main(int argc, char **argv)
{
    auto &app = initialize(argc, argv);
    // app.set_log_level(LogLevel::debug);

    // Mesh, order and FE space
    auto mesh = io::read_mesh(argv[1]);
    const int order = 2;
    const auto V = std::make_shared<fem::CGSpace>(mesh, order);

    // Temperature
    FEField T(V, {"T"});

    // Dirichlet B.C.
    fem::DirichletBC bc(V);
    bc.set_value("Left", 10);
    bc.set_value("Right", 100);

    // Diffusivity
    ConstantField D("D", 1.0);

    auto Axb = fem::create_axb(T, la::SolverType::cg,
                               {}, la::Backend::petsc);
    kernels::Diffusion diffusion(T, D);
    diffusion(Axb->lhs(), Axb->rhs());
    Axb->assemble();

    const auto [dof_idxs, dof_values] = bc.get_dofs_values();
    Axb->eliminate_dofs(dof_idxs, dof_values);
    Axb->solve(T.dof_values());
    T.dof_values().update_ghosts();

    // Write solution to file
    io::vtk::write("post/solution_000", {T});

    return 0;
}