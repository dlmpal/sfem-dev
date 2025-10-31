#pragma once

#include "flux_function.hpp"

namespace sfem::fvm
{
    class EulerFlux : public FluxFunction
    {
    public:
        EulerFlux(real_t gamma, int dim);

        real_t compute_flux(const std::vector<real_t> &state, std::vector<real_t> &flux, int dir) const override;

        real_t compute_normal_flux(const std::vector<real_t> &state,
                                   const geo::Vec3 &normal,
                                   std::vector<real_t> &normal_flux) const override;

    private:
        /// @brief Adiabatic index
        real_t gamma_;
    };
}