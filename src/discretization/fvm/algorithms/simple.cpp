#include "simple.hpp"
#include "../kernels/transient.hpp"
#include "../kernels/convection.hpp"
#include "../kernels/laplacian.hpp"
#include "../kernels/source.hpp"
#include "../utils/la_utils.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fvm::algo
{
    //=============================================================================
    SIMPLESolver::SIMPLESolver(std::vector<FVField> U, FVField P,
                               real_t rho, real_t mu, SIMPLEOptions options)
        : U_(U),
          P_(P),
          Pcorr_(P_.space(), {"Pcorr"}, P_.grad_method()),
          D_(P_.space(), {"D"}),
          rho_("rho", rho),
          mu_("mu", mu),
          pressure_(Pcorr_),
          options_(options),
          dt_(1.0)
    {
        // For the pressure correction field,
        // all Dirichlet BCs are set to 0
        Pcorr_.boundary_condition() = P.boundary_condition();
        const auto mesh = P_.space()->mesh();
        for (const auto &region : mesh->regions())
        {
            if (region.dim() < mesh->pdim())
            {
                if (P_.boundary_condition().region_type(region.name()) == BCType::dirichlet)
                {
                    Pcorr_.boundary_condition().set_region_bc(region.name(), BCType::dirichlet, 0.0);
                }
            }
        }

        /// @todo
        flux_.resize(mesh->topology()->n_entities(mesh->pdim() - 1), 0.0);

        // Setup momentum equation
        auto momentum_Axb = create_axb(U_.front(),
                                       options_.momentum_solver_type,
                                       options_.momentum_solver_options,
                                       options_.backend);
        int dir = 0;
        for (auto u : U_)
        {
            Equation eqn(u, momentum_Axb);

            auto rhs = [P, dir](const FVField &, int cell_idx, std::span<real_t> dpdxi)
            {
                dpdxi[0] = -P.cell_grad(cell_idx)(dir);
            };

            if (options_.transient)
            {
                eqn.add_kernel(ImplicitEuler(u, rho_, dt_));
            }
            eqn.add_kernel(Convection(u, flux_));
            eqn.add_kernel(Laplacian(u, mu_));
            eqn.add_kernel(Source(u, rhs));

            momentum_.emplace_back(std::move(eqn));

            dir++;
        }

        // Setup pressure equation
        auto pressure_Axb = create_axb(Pcorr_,
                                       options_.pressure_solver_type,
                                       options_.pressure_solver_options,
                                       options_.backend);
        pressure_ = Equation(Pcorr_, pressure_Axb);
        auto pressure_rhs = [&](la::MatSet, la::VecSet b)
        {
            const auto V = P_.space();

            auto work = [&](const mesh::Mesh &,
                            const mesh::Region &,
                            const mesh::Cell &,
                            int facet_idx)
            {
                const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
                const auto [owner, neighbour] = adjacent_cells;
                std::array<real_t, 2> values = {};
                if (owner == neighbour)
                {
                    values[0] = -flux_[facet_idx] / rho_.cell_value(facet_idx);
                }
                else
                {
                    values[0] = -flux_[facet_idx] / rho_.facet_value(facet_idx);
                    values[1] = flux_[facet_idx] / rho_.facet_value(facet_idx);
                }
                b(adjacent_cells, values);
            };
            mesh::utils::for_all_facets(*V->mesh(), work);
        };
        pressure_.add_kernel(Laplacian(Pcorr_, D_));
        pressure_.add_kernel(pressure_rhs);
    }
    //=============================================================================
    void SIMPLESolver::step(real_t time, real_t dt)
    {
        dt_ = dt;

        for (auto &eqn : momentum_)
        {
            eqn.assemble();
            eqn.apply_relaxation(options_.momentum_alpha);
            eqn.solve();
        }

        assemble_pressure_diffusivity();
        assemble_mass_flux();

        Pcorr_.values().set_all(0.0);
        for (int iter = 0; iter < options_.n_orthogonal_correctors + 1; iter++)
        {
            pressure_.assemble();
            pressure_.solve();
            Pcorr_.update_gradient();
        }

        correct_fields();
    }
    //=============================================================================
    void SIMPLESolver::assemble_pressure_diffusivity()
    {
        const auto V = P_.space();
        const auto mesh = V->mesh();
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &,
                        int cell_idx)
        {
            real_t a_avg = 0.0;
            for (int dir = 0; dir < mesh->pdim(); dir++)
            {
                a_avg += momentum_[dir].diag()(cell_idx);
            }
            a_avg = a_avg / mesh->pdim();
            D_.cell_value(cell_idx) = V->cell_volume(cell_idx) / a_avg;
        };
        mesh::utils::for_all_cells(*mesh, work);
        D_.values().update_ghosts();
    }
    //=============================================================================
    void SIMPLESolver::assemble_mass_flux()
    {
        const auto V = P_.space();
        const auto mesh = V->mesh();
        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &region,
                        const mesh::Cell &,
                        int facet_idx)
        {
            const auto [owner, neighbour] = V->facet_adjacent_cells(facet_idx);

            const real_t Af = V->facet_area(facet_idx);
            const geo::Vec3 nf = V->facet_normal(facet_idx);
            const real_t d12 = V->intercell_distance(facet_idx).mag();

            real_t Uf = 0.0;
            real_t rhof = 0.0;
            real_t Df = 0.0;
            real_t gradPf = 0.0;
            if (owner != neighbour)
            {
                for (int dir = 0; dir < mesh->pdim(); dir++)
                {
                    Uf += U_[dir].facet_value(facet_idx) * nf(dir);
                }
                rhof = rho_.facet_value(facet_idx);
                Df = D_.facet_value(facet_idx);
                /// @todo Implement Rhie-Chow correctly
                gradPf = (P_.cell_value(neighbour) - P_.cell_value(owner)) / d12 - geo::inner(P_.facet_grad(facet_idx), nf);
                flux_[facet_idx] = rhof * (Uf - Df * gradPf) * Af;
            }
            else
            {
                for (int dir = 0; dir < mesh->pdim(); dir++)
                {
                    const FVBC &ubc = U_[dir].boundary_condition();
                    real_t uf = U_[dir].cell_value(owner);
                    if (ubc.region_type(region.name()) == fvm::BCType::dirichlet)
                    {
                        uf = ubc.facet_value(facet_idx);
                    }
                    Uf += uf * nf(dir);
                }
                rhof = rho_.cell_value(owner);
                flux_[facet_idx] = rhof * Uf * Af;
            }
        };
        mesh::utils::for_all_facets(*mesh, work);
    }
    //=============================================================================
    void SIMPLESolver::correct_fields()
    {
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
            mesh::utils::for_all_cells(*mesh, work, false);

            for (auto &u : U_)
            {
                u.update_gradient();
            }
        }

        // Correct mass fluxes
        {
            auto work = [&](const mesh::Mesh &,
                            const mesh::Region &,
                            const mesh::Cell &,
                            int facet_idx)
            {
                const real_t Af = V->facet_area(facet_idx);
                const geo::Vec3 nf = V->facet_normal(facet_idx);

                const real_t rhof = rho_.facet_value(facet_idx);
                const real_t Df = D_.facet_value(facet_idx);
                const geo::Vec3 gradPf = Pcorr_.facet_grad(facet_idx);

                flux_[facet_idx] -= rhof * Df * Af * geo::inner(gradPf, nf);
            };
            mesh::utils::for_all_facets(*mesh, work);
        }

        // Correct pressure and compute gradient
        la::axpy(options_.pressure_alpha, Pcorr_.values(), P_.values());
        P_.values().update_ghosts();
        P_.update_gradient();
    }
}