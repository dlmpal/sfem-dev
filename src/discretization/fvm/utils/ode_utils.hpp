#pragma once

#include "../fv_field.hpp"
#include "../hyperbolic.hpp"
#include "../../ode/erk.hpp"

namespace sfem::fvm::ode
{
    using RHSFunction = sfem::ode::ERKIntegrator::RHSFunction;

    /// @brief Create a RHS function for a given finite volume field, numerical flux and source function
    RHSFunction create_rhs(std::shared_ptr<const fvm::FVField> phi,
                           std::shared_ptr<const fvm::NumericalFlux> nflux,
                           FieldFunction src = {});
}