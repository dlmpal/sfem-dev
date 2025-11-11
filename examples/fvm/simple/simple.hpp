#include "sfem.hpp"

namespace sfem::fvm::simple
{
    /// @brief Assemble the momentum equation LHS matrix and RHS vector for a given direction
    /// @param U Velocity field for the given direction
    /// @param gradU Velocity field gradient
    /// @param bc Boundary conditions
    /// @param gradP Pressure field gradient
    /// @param nu  Kinematic viscosity
    /// @param vel Convective velocity
    /// @param rho Density
    /// @param dt Timestep size
    /// @param dir Direction
    /// @param A Matset for LHS matrix
    /// @param b Vecset for RHS vector
    void assemble_momentum(const FVField &U,
                           const FVField &gradU,
                           const FVBC &bc,
                           const FVField &gradP,
                           const Coefficient &nu,
                           const Coefficient &vel,
                           real_t rho, real_t dt, int dir,
                           la::MatSet A, la::VecSet b)
    {
        laplacian(U, gradU, bc, nu, A, b);
        convection(U, bc, vel, A, b);
        if (dt > 0)
        {
            ddt(U, ConstantCoefficient(1.0), dt, A, b);
        }

        // Add pressure gradient to RHS
        const auto V = U.space();
        const real_t rho_inv = 1.0 / rho;
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            std::array<int, 1> idx = {cell_idx};
            std::array<real_t, 1> value = {-rho_inv * gradP(cell_idx, dir) *
                                           V->cell_volume(cell_idx)};
            b(idx, value);
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }

    /// @brief Assemble the pressure equation LHS matrix and RHS vector
    /// @param P Pressure field
    /// @param gradP Pressure field gradient
    /// @param bc Pressure boundary conditions
    /// @param D Diffusivity
    /// @param a Momentum matrix diagonal coefficiets
    /// @param H Momentum matrix off-diagonal coefficients
    /// @param A Matset for LHS matrix
    /// @param b Vecset for RHS vector
    void assemble_pressure(const FVField &P,
                           const FVField &gradP,
                           const FVBC &bc,
                           const Coefficient &D,
                           const la::Vector &a,
                           const la::Vector &H,
                           la::MatSet A, la::VecSet b)
    {
        // Quick access
        const auto V = P.space();
        const int dim = V->mesh()->pdim();

        // Assemble Laplacian
        laplacian(P, gradP, bc, D, A, b);

        // Assemble RHS
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int facet_idx)
        {
            const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
            const auto [cell_idx1, cell_idx2] = adjacent_cells;

            const real_t g = V->facet_interp_factor(facet_idx);
            const real_t Af = V->facet_area(facet_idx);
            const geo::Vec3 nf = V->facet_normal(facet_idx);

            std::array<real_t, 2> values = {};
            for (int i = 0; i < dim; i++)
            {
                const real_t a1 = a(cell_idx1, i);
                const real_t a2 = a(cell_idx2, i);

                const real_t H1 = H(cell_idx1, i);
                const real_t H2 = H(cell_idx2, i);

                if (cell_idx1 != cell_idx2)
                {
                    const real_t avg = g * (H1 / a1) + (1 - g) * (H2 / a2);
                    values[0] -= avg * nf(i) * Af;
                    values[1] += avg * nf(i) * Af;
                }
                else
                {
                    values[0] -= H1 / a1 * nf(i) * Af;
                    values[1] = 0.0;
                }
            }
            b(adjacent_cells, values);
        };
        mesh::utils::for_all_facets(*V->mesh(), work);
    }

    /// @brief Apply explicit relaxation to field
    /// @param phi_old Old field
    /// @param phi_new New field (to be relaxed)
    /// @param factor Relaxation factor
    void apply_relaxation(const FVField &phi_old,
                          FVField &phi_new,
                          real_t factor)
    {
        la::scale(factor, phi_new);
        la::axpy(1 - factor, phi_old, phi_new);
    }

    /// @brief Compute the mass residual based on the continuity equation: div(U)=0
    /// @param gradU Velocity field gradient
    /// @return Mass residual
    real_t compute_mass_residual(const std::vector<std::shared_ptr<FVField>> &gradU)
    {
        const auto V = gradU[0]->space();
        const int dim = V->mesh()->pdim();

        real_t mass_residual = 0.0;
        real_t total_volume = 0.0;
        for (int i = 0; i < V->index_map()->n_owned(); ++i)
        {
            real_t cell_div = 0.0;
            for (int j = 0; j < dim; j++)
            {
                cell_div += (*gradU[j])(i);
            }
            cell_div *= V->cell_volume(i);

            mass_residual += std::abs(cell_div);
            total_volume += V->cell_volume(i);
        }

        mass_residual = mpi::reduce(mass_residual, mpi::ReduceOperation::sum);
        total_volume = mpi::reduce(total_volume, mpi::ReduceOperation::sum);

        return mass_residual / total_volume;
    }
}