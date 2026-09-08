#pragma once

#include <sfem/discretization/fvm/core/fv_equation.hpp>

namespace sfem::fvm::simple
{
    void setup_momentum_eqn(FVField u, FVField P_, Field &mu,
                            std::vector<real_t> &mdot,
                            Equation &eqn, int dir);

    void setup_pressure_eqn(FVField &P, Field &D, Field &rho,
                            const std::vector<real_t> &mdot,
                            Equation &eqn);

    void compute_mass_flux(const std::vector<FVField> &U,
                           const FVField &P,
                           const FVField &D,
                           const Field &rho,
                           std::vector<real_t> &mdot);

    real_t compute_mass_residual(std::shared_ptr<const FVSpace> V,
                                 const std::vector<real_t> &mdot);
}