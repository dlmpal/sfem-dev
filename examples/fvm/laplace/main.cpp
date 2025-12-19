// Solve Laplace's equation on a 2D or 3D domain
// Create the mesh with the cart-mesh application by:
// ${SFEM_DEV_INSTALL_DIR}/bin/cart-mesh -d=2 -Nx=100 -Ny=100 -x-low=0 -x-high=1 -y-low=0 -y-high=1

#include <sfem/sfem.hpp>

using namespace sfem;
using namespace fvm;

int main(int argc, char **argv)
{
    initialize(argc, argv);

    // Read the mesh from file
    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);

    // Finite volume space
    auto V = std::make_shared<FVSpace>(mesh);

    // Temperature field
    FVField T(V, {"T"}, GradientMethod::green_gauss);

    // Boundary conditions
    T.boundary_condition().set_region_bc("Left", BCType::dirichlet, 1.0);
    T.boundary_condition().set_region_bc("Right", BCType::dirichlet, 0.0);

    // Diffusivity
    ConstantField D("D", 1.0);

    // Equation
    Equation eqn(T);
    eqn.add_kernel(Laplacian(T, D));

    const int n_orthogonal_correctors = 1;
    for (int i = 0; i < n_orthogonal_correctors + 1; i++)
    {
        log_msg(std::format("Non-Orthogonal Corrector - Iteration: {}\n", i), true);

        eqn.assemble();
        eqn.solve();

        // Save solution to file
        io::vtk::write(std::format("post/solution_00{}", i), {T});
    }

    return 0;
}