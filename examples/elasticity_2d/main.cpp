#include "sfem.hpp"
#include <iostream>

using namespace sfem;

class LinearElasticity2D
{
public:
    LinearElasticity2D(real_t E, real_t nu, real_t thick)
        : E_(E), nu_(nu), thick_(thick)
    {
    }

    la::DenseMatrix operator()(int cell_idx, const fem::FEData &data)
    {
        // Stress-strain matrix
        la::DenseMatrix D(3, 3);
        real_t coeff = thick_ * E_ / (1 - nu_ * nu_);
        D(0, 0) = coeff * 1.0;
        D(0, 1) = coeff * nu_;
        D(1, 0) = coeff * nu_;
        D(1, 1) = coeff * 1.0;
        D(2, 2) = coeff * (1 - nu_) * 0.5;

        // Strain-displacement matrix
        const auto &dNdX = data.dNdX;
        int n_cols = dNdX.n_rows() * 2;
        int n_rows = 3;
        la::DenseMatrix B(n_rows, n_cols);
        for (int i = 0; i < dNdX.n_rows(); i++)
        {
            B(0, i * 2 + 0) = dNdX(i, 0); ///< exx
            B(1, i * 2 + 1) = dNdX(i, 1); ///< eyy
            B(2, i * 2 + 0) = dNdX(i, 1); ///< exy
            B(2, i * 2 + 1) = dNdX(i, 0); ///< exy
        }

        return B.T() * D * B;
    }

private:
    real_t E_;
    real_t nu_;
    real_t thick_;
};

class PressureLoad2D
{
public:
    PressureLoad2D(real_t thick, real_t pressure_value)
        : thick_(thick), pressure_value_(pressure_value)
    {
    }

    la::DenseMatrix operator()(int facet_idx,
                               const fem::FEData &data,
                               const geo::Vec3 &normal)
    {
        la::DenseMatrix F(data.N.n_rows() * 2, 1, 0.0);
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            F(i * 2 + 0, 0) += -pressure_value_ * thick_ * data.N(i, 0) * normal.x();
            F(i * 2 + 1, 0) += -pressure_value_ * thick_ * data.N(i, 0) * normal.y();
        }
        return F;
    }

private:
    real_t thick_;
    real_t pressure_value_;
};

int main(int argc, char **argv)
{
    initialize(argc, argv);

    auto mesh = io::read_mesh(argv[1]);
    int dim = mesh->topology()->dim();

    const int order = 1;
    auto phi = std::make_shared<fem::CGSpace>(mesh, order, std::vector<std::string>{"Ux", "Uy"});

    const real_t thick = 0.5 * 1e-3;
    const real_t E = 5e9;
    const real_t nu = 0.35;
    const real_t pressure_value = 1000;

    // LHS matrix
    auto A = fem::petsc::create_mat(*phi);
    fem::assemble_matrix_cells(*phi, "Solid", LinearElasticity2D(E, nu, thick), fem::petsc::create_matset(A));
    A.assemble();

    // RHS vector
    auto b = fem::petsc::create_vec(*phi);
    fem::assemble_vec_facets(*phi, "Left", PressureLoad2D(thick, pressure_value), fem::petsc::create_vecset(b));
    b.assemble();

    // Dirichlet B.C.
    fem::DirichletBC bc(phi);
    bc.set_value("Fixed", "Ux", 0);
    bc.set_value("Fixed", "Uy", 0);

    // Solution vector
    auto x = fem::petsc::create_vec(*phi);

    // Apply Dirichlet BC and solve
    fem::petsc::solve(A, b, x, bc);

    // Save solution to VTK file
    io::vtk::write("post/solution_000", *phi, x.get_values());

    return 0;
}