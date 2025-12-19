#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <sfem/discretization/fvm/physics/hyperbolic/numerical_flux.hpp>
#include <sfem/discretization/ode/erk.hpp>

namespace sfem::fvm::ode
{
    using RHSFunction = sfem::ode::ERKIntegrator::RHSFunction;

    /// @brief Create a RHS function for a given finite volume field, numerical flux and source function
    RHSFunction create_rhs(const FVField &phi,
                           std::shared_ptr<const fvm::NumericalFlux> nflux);
}