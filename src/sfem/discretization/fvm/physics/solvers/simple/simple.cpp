#include "simple.hpp"
#include <sfem/discretization/fvm/physics/solvers/simple/utils.hpp>
#include <sfem/discretization/fvm/physics/kernels/transient.hpp>
#include <sfem/discretization/fvm/physics/kernels/convection.hpp>
#include <sfem/discretization/fvm/physics/kernels/laplacian.hpp>
#include <sfem/discretization/fvm/physics/kernels/source.hpp>
#include <sfem/discretization/fvm/core/utils/la_utils.hpp>
#include <sfem/mesh/utils/loop_utils.hpp>

namespace sfem::fvm::simple
{
    //=============================================================================
    SimpleSolver::SimpleSolver(const std::vector<FVField> &U, const FVField &P,
                               real_t rho, real_t mu, SimpleOptions options)
        : U_(U),
          P_(P),
          Pcorr_(P.space(), {"Pcorr"}, P.grad_method()),
          D_(P.space(), {"D"}),
          rho_("rho", rho),
          mu_("mu", mu),
          pressure_(Pcorr_),
          options_(options)
    {
        // For the pressure correction field Dirichlet BCs are set to 0
        Pcorr_.boundary_condition() = P.boundary_condition();
        const auto mesh = P_.space()->mesh();
        for (const auto &region : mesh->regions())
        {
            if (region.dim() < mesh->pdim())
            {
                if (Pcorr_.boundary_condition().region_type(region.name()) == BCType::dirichlet)
                {
                    Pcorr_.boundary_condition().set_region_bc(region.name(), BCType::dirichlet, 0.0);
                }
            }
        }

        // Allocate mass flux vector
        mdot_.resize(mesh->topology()->n_entities(mesh->pdim() - 1), 0.0);

        // Setup momentum equation
        auto momentum_Axb = create_axb(U_.front(),
                                       options_.momentum_solver_type,
                                       options_.momentum_solver_options,
                                       options_.backend);
        int dir = 0;
        for (auto u: U_)
        {
            Equation eqn(u, momentum_Axb);
            if (options_.transient)
            {
                U_old_.push_back(FVField(u.space(), u.components()));
                eqn.add_kernel(ImplicitEuler(U_old_[dir], rho_, dt_));
            }
            setup_momentum_eqn(u, P_, mu_, mdot_, eqn, dir++);
            momentum_.emplace_back(std::move(eqn));
        }

        // Setup pressure equation
        auto pressure_Axb = create_axb(Pcorr_,
                                       options_.pressure_solver_type,
                                       options_.pressure_solver_options,
                                       options_.backend);
        pressure_ = Equation(Pcorr_, pressure_Axb);
        setup_pressure_eqn(Pcorr_, D_, rho_, mdot_, pressure_);

        plot_file_ = options_.plot_file;
        plot_int_ = options.plot_int;
        plot_fields_ = U_;
        plot_fields_.push_back(P_);
    }
    //=============================================================================
    void SimpleSolver::pre_timestep()
    {
        if (options_.transient)
        {
            for (std::size_t i = 0; i < U_.size(); i++)
            {
                la::copy(U_[i].values(), U_old_[i].values());
                U_old_[i].values().update_ghosts();
            }
        }
    }
    //=============================================================================
    void SimpleSolver::solve_momentum()
    {
        for (auto& eqn: momentum_)
        {
            eqn.assemble();
            eqn.apply_relaxation(options_.momentum_alpha);
            eqn.solve();
            residual_history_.push_back(eqn.Axb()->residual_history().front());
        }
    }
    //=============================================================================
    void SimpleSolver::solve_pressure()
    {
        // Update diffusivity
        const auto V = P_.space();
        const auto mesh = V->mesh();
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            real_t a_avg = 0.0;
            for (const auto& eqn: momentum_)
            {
                a_avg += eqn.diag()(cell_idx);
            }
            a_avg = a_avg / mesh->pdim();
            D_.cell_value(cell_idx) = V->cell_volume(cell_idx) / a_avg;
        };
        mesh::utils::for_all_cells(*mesh, work);
        D_.values().update_ghosts();

        // Update mass flux
        compute_mass_flux(U_, P_, D_, rho_, mdot_);

        Pcorr_.values().set_all(0.0);
        Pcorr_.update_gradient();
        for (int iter = 0; iter < options_.n_orthogonal_correctors + 1; iter++)
        {
            pressure_.assemble();
            pressure_.solve();
        }
        residual_history_.push_back(pressure_.Axb()->residual_history().front());
    }
    //=============================================================================
    void SimpleSolver::correct_fields()
    {
        // Quick access
        const auto V = P_.space();
        const auto mesh = V->mesh();

        // Correct velocity
        {
            auto work = [&](const mesh::Mesh &,
                            const mesh::Region &,
                            const mesh::Cell &,
                            int cell_idx)
            {
                for (int dir = 0; dir < mesh->pdim(); dir++)
                {
                    U_[dir].cell_value(cell_idx) -= D_.cell_value(cell_idx) * Pcorr_.cell_grad(cell_idx)(dir);
                }
            };
            mesh::utils::for_all_cells(*mesh, work);

            for (auto &u : U_)
            {
                u.values().update_ghosts();
                u.update_gradient();
            }
        }

        // Correct mass flux
        {
            const FVBC &pcorr_bc = Pcorr_.boundary_condition();

            auto work = [&](const mesh::Mesh &,
                            const mesh::Region &region,
                            const mesh::Cell &,
                            int facet_idx)
            {
                const auto [owner, neighbour] = V->facet_adjacent_cells(facet_idx);
                const geo::Vec3 dPN = V->facet_intercell_distance(facet_idx);

                // Boundary facets
                if (owner == neighbour)
                {
                    // The flux can only be corrected where the pressure
                    // correction is prescribed. Everywhere else (walls, inlets)
                    // the flux is imposed by the velocity B.C. and the pressure
                    // equation contributes nothing to this facet
                    if (pcorr_bc.region_type(region.name()) != BCType::dirichlet)
                    {
                        return;
                    }

                    const geo::Vec3 Sf = V->facet_area_vec(facet_idx);
                    const real_t rhof = rho_.cell_value(owner);
                    const real_t Df = D_.cell_value(owner);
                    mdot_[facet_idx] += rhof * Df * 2.0 * Sf.mag() / dPN.mag() *
                                        (Pcorr_.cell_value(owner) - pcorr_bc.value(facet_idx));
                }
                // Internal facets
                else
                {
                    const auto [delta, kappa] = V->decompose_area_vec(facet_idx);
                    const real_t rhof = rho_.facet_value(facet_idx);
                    const real_t Df = D_.facet_value(facet_idx);

                    // Orthogonal contribution
                    mdot_[facet_idx] += rhof * Df * delta.mag() / dPN.mag() *
                                        (Pcorr_.cell_value(owner) - Pcorr_.cell_value(neighbour));

                    // Non-orthogonal correction
                    mdot_[facet_idx] -= rhof * Df * geo::inner(Pcorr_.facet_grad(facet_idx), kappa);
                }
            };
            mesh::utils::for_all_facets(*mesh, work);
        }

        // Correct pressure and compute gradient
        la::axpy(options_.pressure_alpha, Pcorr_.values(), P_.values());
        P_.values().update_ghosts();
        P_.update_gradient();
    }
    //=============================================================================
    void SimpleSolver::step(real_t time, real_t dt)
    {
        time_ = time;
        dt_ = dt;

        // Initial mass residual
        real_t mass_residual_0 = 0.0;

        for (int iter = 0; iter < options_.max_iter_simple; iter++)
        {
            solve_momentum();

            solve_pressure();

            correct_fields();

            const real_t mass_residual = compute_mass_residual(P_.space(), mdot_);

            if (iter == 0)
            {
                mass_residual_0 = mass_residual;
            }

            const real_t rel_residual = mass_residual / mass_residual_0;
            log_msg(std::format("SIMPLE Iteration: {}, mass residual: {:.6e}\n",
                                iter, mass_residual),
                    true);

            // Check for convergence based on the mass residual
            if (mass_residual < options_.atol_simple || rel_residual < options_.rtol_simple)
            {
                log_msg(std::format("SIMPLE converged in {} iterations\n", iter + 1), true);
                break;
            }
        }
    }
}