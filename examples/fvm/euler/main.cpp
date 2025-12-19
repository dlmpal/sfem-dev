#include "sfem.hpp"

using namespace sfem;

void riemann1d(fvm::FVField &U, real_t gamma);
void riemann2d1(fvm::FVField &U, real_t gamma);
void riemann2d2(fvm::FVField &U, real_t gamma);
void cylindrical(fvm::FVField &U, real_t gamma);

namespace sfem::fvm::ode
{
    double compute_dt_estimate(const fvm::FVField &U,
                               const fvm::FluxFunction &flux,
                               double CFL)
    {
        SFEM_CHECK_SIZES(U.n_comp(), flux.n_comp());

        // Quick access
        const int dim = flux.dim();
        const int n_comp = flux.n_comp();
        const auto V = U.space();

        // Cell state and flux
        std::vector<real_t> u1(n_comp);
        std::vector<real_t> u2(n_comp);
        std::vector<real_t> f1(n_comp);
        std::vector<real_t> f2(n_comp);

        real_t dt = std::numeric_limits<real_t>::max();

        auto work = [&](const mesh::Mesh &mesh,
                        const mesh::Region &region,
                        const mesh::Cell &facet, int facet_idx)
        {
            const auto [cell_idx1, cell_idx2] = V->facet_adjacent_cells(facet_idx);
            const geo::Vec3 normal = V->facet_normal(facet_idx);
            const real_t d12 = V->intercell_distance(facet_idx).mag();

            for (int n = 0; n < n_comp; n++)
            {
                u1[n] = U(cell_idx1, n);
                u2[n] = U(cell_idx2, n);
            }

            const real_t s1 = flux.compute_normal_flux(u1, normal, f1);
            const real_t s2 = flux.compute_normal_flux(u2, normal, f2);
            const real_t s = std::max(s1, s2);

            const real_t dt_facet = CFL * d12 / s;
            if (dt_facet < dt)
            {
                dt = dt_facet;
            }
        };
        mesh::utils::for_all_facets(*V->mesh(), work, true, true);

        return mpi::reduce(dt, mpi::ReduceOperation::min);
    }

    double minbee(double r)
    {
        if (r <= 0)
        {
            return 0;
        }
        else if (r < 1)
        {
            return r;
        }
        else
        {
            return std::min(1.0, 2 / (1 + r));
        }
    }

    double van_leer(double r)
    {
        return (r + std::abs(r)) / (1 + std::abs(r));
    }

    RHSFunction create_rhs_(std::shared_ptr<const FVField> phi,
                            std::shared_ptr<const FVField> grad,
                            std::shared_ptr<const NumericalFlux> nflux,
                            FieldFunction src = {})
    {
        return [=](const la::Vector &S,
                   la::Vector &rhs,
                   real_t time)
        {
            // Reset RHS vector
            rhs.set_all(0.0);

            // Finite volume space and flux function
            const auto V = phi->space();
            const auto flux = nflux->flux_function();

            fvm::FVBC bc(phi);

            // Left and right state fluxes and normal flux
            std::vector<real_t> state_left(flux->n_comp());
            std::vector<real_t> state_right(flux->n_comp());
            std::vector<real_t> normal_flux(flux->n_comp());

            // Add flux contribution
            auto facet_work = [&](const mesh::Mesh &,
                                  const mesh::Region &,
                                  const mesh::Cell &,
                                  int facet_idx)
            {
                // Get the facet's adjacent cells
                const auto [cell_idx_left, cell_idx_right] = V->facet_adjacent_cells(facet_idx);

                // Facet area and normal vector
                const real_t area = V->facet_area(facet_idx);
                const auto normal = V->facet_normal(facet_idx);
                const real_t g = V->facet_interp_factor(facet_idx);
                const geo::Vec3 d = V->intercell_distance(facet_idx).normalize();

                // Get the adjacent cell states
                for (int i = 0; i < flux->n_comp(); i++)
                {
                    state_left[i] = S(cell_idx_left, i);
                    state_right[i] = S(cell_idx_right, i);
                }

                // Reconstruct the adjacent cell states
                real_t normal_vel = 0.0;
                for (int i = 0; i < flux->dim(); i++)
                {
                    normal_vel += (g * state_left[i + 1] + (1 - g) * state_right[i + 2]) * normal(i);
                }

                const int upwind_idx = normal_vel > 0 ? cell_idx_left : cell_idx_right;
                const int downwind_idx = normal_vel < 0 ? cell_idx_left : cell_idx_right;

                const geo::Vec3 x_left(V->cell_midpoint(cell_idx_left),
                                       V->facet_midpoint(facet_idx));
                const geo::Vec3 x_right(V->cell_midpoint(cell_idx_right),
                                        V->facet_midpoint(facet_idx));

                const real_t eps = 1e-10;
                for (int i = 0; i < flux->n_comp(); i++)
                {
                    real_t denom = 0.0;
                    for (int j = 0; j < flux->dim(); j++)
                    {
                        denom += (*grad)(upwind_idx, i * flux->dim() + j) * d(j);
                    }
                    //  const real_t denom = geo::inner(gr)
                    const real_t r = 1.0 - 0.5 * (S(downwind_idx, i) - S(upwind_idx, i)) / (denom + eps);

                    for (int j = 0; j < flux->dim(); j++)
                    {
                        state_left[i] += 0.5 * van_leer(r) * x_left(j) * (*grad)(cell_idx_left, i * flux->dim() + j);
                        state_right[i] += 0.5 * van_leer(r) * x_right(j) * (*grad)(cell_idx_right, i * flux->dim() + j);
                    }
                }

                // Compute the (numerical) normal flux at the face
                nflux->compute_normal_flux(state_left, state_right, normal, normal_flux);

                // Inverse left and right cell volumes
                const real_t vol_inv_left = 1 / V->cell_volume(cell_idx_left);
                const real_t vol_inv_right = 1 / V->cell_volume(cell_idx_right);

                // Add the flux contributions to the RHS vector
                for (int i = 0; i < flux->n_comp(); i++)
                {
                    rhs(cell_idx_left, i) -= normal_flux[i] * area * vol_inv_left;

                    /// @todo cleanup
                    if (cell_idx_left != cell_idx_right)
                    {
                        rhs(cell_idx_right, i) -= -normal_flux[i] * area * vol_inv_right;
                    }
                }
            };
            mesh::utils::for_all_facets(*V->mesh(), facet_work);

            // Add source term contribution
            if (src)
            {
                std::vector<real_t> source(flux->n_comp());
                auto cell_work = [&](const mesh::Mesh &,
                                     const mesh::Region &,
                                     const mesh::Cell &,
                                     int cell_idx)
                {
                    src(V->cell_midpoint(cell_idx), source, time);
                    for (int i = 0; i < flux->n_comp(); i++)
                    {
                        rhs(cell_idx, i) += source[i];
                    }
                };
                mesh::utils::for_all_cells(*V->mesh(), cell_work);
            };

            // RHS was incrementally constructed - needs assembly
            rhs.assemble();
        };
    }
}

int main(int argc, char *argv[])
{
    initialize(argc, argv);

    // Read the mesh from file
    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);

    // Finite volume space
    auto V = std::make_shared<fvm::FVSpace>(mesh);

    // State vectors, gradient and RHS
    auto U_new = fvm::create_field(V, {"rho", "momx", "momy", "E"});
    auto U_old = fvm::create_field(V, {"rho", "momx", "momy", "E"});
    auto gradU = fvm::create_gradient(*U_new);
    auto rhs = fvm::create_field(V, {"rho", "momx", "momy", "E"});

    // Adiabatic index
    const real_t gamma = 1.4;

    // Initial condition
    cylindrical(*U_old, gamma);
    U_old->update_ghosts();
    io::vtk::write(std::format("post/solution_{}", 0), *mesh, {U_old});

    // Flux function and numerical flux
    auto flux = std::make_shared<fvm::EulerFlux>(gamma, 2);
    auto nflux = std::make_shared<fvm::RusanovFlux>(flux);

    // Integrator
    auto rhs_func = fvm::ode::create_rhs_(U_new, gradU, nflux);
    auto erk = ode::create_erk(*U_new, rhs_func, ode::ERKType::rk2);

    const real_t CFL = 0.35;
    real_t dt = fvm::ode::compute_dt_estimate(*U_old, *flux, CFL);
    const real_t time_start = 0.;
    const real_t time_stop = 0.25;
    const int plot_int = 10;
    int step = 0;
    real_t time = time_start;
    while (time <= time_stop)
    {
        // Increment time
        step++;
        time += dt;
        log_msg(std::format("Time: {}, Step: {}\n", time, step), true);

        // Advance
        erk.advance(*U_old, *U_new, time, dt);

        // Recompute timestep
        real_t dt_ = fvm::ode::compute_dt_estimate(*U_new, *flux, CFL);
        if (dt_ < dt)
        {
            dt = dt_;
        }

        // Save solution to file
        if (step % plot_int == 0)
        {
            U_new->update_ghosts();
            io::vtk::write(std::format("post/solution_{}", step), *mesh, {U_new});
        }

        la::copy(*U_new, *U_old);
    }

    return 0;
}

void riemann1d(fvm::FVField &U, real_t gamma)
{
    const int dim = U.n_comp() - 2;
    const auto V = U.space();

    const double rho_l = 1.0;
    const double u_l = 0.0;
    const double p_l = 1.0;
    const double e_l = p_l / rho_l / (gamma - 1);
    const double E_l = rho_l * (e_l + 0.5 * u_l * u_l);

    const double rho_r = 0.125;
    const double u_r = 0.0;
    const double p_r = 0.1;
    const double e_r = p_r / rho_r / (gamma - 1);
    const double E_r = rho_r * (e_r + 0.5 * u_r * u_r);

    auto work = [&](const mesh::Mesh &mesh,
                    const mesh::Region &region,
                    const mesh::Cell &cell, int cell_idx)
    {
        if (V->cell_midpoint(cell_idx)[1] < 0.5)
        {
            U(cell_idx, 0) = rho_l;
            U(cell_idx, 2) = rho_l * u_l;
            U(cell_idx, dim + 1) = E_l;
        }
        else
        {
            U(cell_idx, 0) = rho_r;
            U(cell_idx, 2) = rho_r * u_r;
            U(cell_idx, dim + 1) = E_r;
        }
    };
    mesh::utils::for_all_cells(*V->mesh(), work);
}

void riemann2d1(fvm::FVField &U, real_t gamma)
{
    const int dim = U.n_comp() - 2;
    const auto V = U.space();

    const real_t rho_high = 1.0;
    const real_t rho_medium = 0.4;
    const real_t rho_low = 0.125;

    const real_t p_high = 1.0;
    const real_t p_medium = 0.3;
    const real_t p_low = 0.1;

    const real_t x_lim = 0.75;
    const real_t y_lim = 0.75;

    auto work = [&](const mesh::Mesh &mesh,
                    const mesh::Region &region,
                    const mesh::Cell &cell,
                    int cell_idx)
    {
        const auto [x, y, z] = V->cell_midpoint(cell_idx);

        if (x > x_lim)
        {
            if (y > y_lim)
            {
                U(cell_idx, 0) = rho_high;
                U(cell_idx, dim + 1) = p_high / (gamma - 1.0);
            }
            else
            {
                U(cell_idx, 0) = rho_medium;
                U(cell_idx, dim + 1) = p_medium / (gamma - 1.0);
            }
        }
        else
        {
            if (y > y_lim)
            {
                U(cell_idx, 0) = rho_medium;
                U(cell_idx, dim + 1) = p_medium / (gamma - 1.0);
            }
            else
            {
                U(cell_idx, 0) = rho_low;
                U(cell_idx, dim + 1) = p_low / (gamma - 1.0);
            }
        }
    };
    mesh::utils::for_all_cells(*V->mesh(), work);
}

void riemann2d2(fvm::FVField &U, real_t gamma)
{
    const int dim = U.n_comp() - 2;
    const auto V = U.space();
    const real_t x_lim = 1.0;
    const real_t y_lim = 1.0;

    auto work = [&](const mesh::Mesh &mesh,
                    const mesh::Region &region,
                    const mesh::Cell &cell,
                    int idx)
    {
        const auto [x, y, z] = V->cell_midpoint(idx);

        if (x > x_lim)
        {
            if (y > y_lim)
            {
                U(idx, 0) = 1.5;
                U(idx, dim + 1) = 1.5 / (gamma - 1.0);
            }
            else
            {
                U(idx, 0) = 33. / 62.;
                U(idx, 2) = 4. / std::sqrt(11);
                U(idx, dim + 1) = 0.3 / (gamma - 1.0);
            }
        }
        else
        {
            if (y > y_lim)
            {
                U(idx, 0) = 33. / 62.;
                U(idx, 1) = 4 / std::sqrt(11);
                U(idx, dim + 1) = 0.3 / (gamma - 1.0);
            }
            else
            {
                U(idx, 0) = 77. / 558.;
                U(idx, 1) = 4 / std::sqrt(11);
                U(idx, 2) = 4 / std::sqrt(11);
                U(idx, dim + 1) = 9. / 310. / (gamma - 1.0);
            }
        }
        //
    };
}

void cylindrical(fvm::FVField &U, real_t gamma)
{
    const int dim = U.n_comp() - 2;
    const auto V = U.space();

    const double rho_l = 1.0;
    const double u_l = 0.0;
    const double v_l = 0.0;
    const double p_l = 1.0;
    const double e_l = p_l / rho_l / (gamma - 1);
    const double E_l = rho_l * e_l;

    const double rho_r = 0.125;
    const double u_r = 0.0;
    const double v_r = 0.0;
    const double p_r = 0.1;
    const double e_r = p_r / rho_r / (gamma - 1);
    const double E_r = rho_r * e_r;

    auto work = [&](const mesh::Mesh &mesh,
                    const mesh::Region &region,
                    const mesh::Cell &cell, int cell_idx)
    {
        const real_t r = geo::compute_distance({}, V->cell_midpoint(cell_idx));
        if (r < 0.4)
        {
            U(cell_idx, 0) = rho_l;
            U(cell_idx, dim + 1) = E_l;
        }
        else
        {
            U(cell_idx, 0) = rho_r;
            U(cell_idx, dim + 1) = E_r;
        }
    };
    mesh::utils::for_all_cells(*V->mesh(), work);
}
