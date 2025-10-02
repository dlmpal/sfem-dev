#pragma once

#include "../fv_function.hpp"
#include "../hyperbolic.hpp"
#include "../../ode/erk.hpp"

namespace sfem::fvm::ode
{
    using SourceFunction = std::function<void(const std::array<real_t, 3> &, std::vector<real_t> &, real_t)>;
    using RHSFunction = sfem::ode::ERKIntegrator::RHSFunction;

    /// @brief Create a RHS function for a given FV function, numerical flux and source function
    RHSFunction create_rhs(std::shared_ptr<const fvm::FVFunction> phi,
                           std::shared_ptr<const fvm::NumericalFlux> nflux,
                           SourceFunction src);
}