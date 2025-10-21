#pragma once

#include "flux_function.hpp"

namespace sfem::fvm
{
    /// @brief Flux function for the two-dimensional Euler equations
    /// @note Subject to change
    class Euler2D : public FluxFunction
    {
    public:
        Euler2D(real_t gamma);

        real_t compute_flux(const std::vector<real_t> &state, std::vector<real_t> &flux, int dir) const override;

        real_t compute_normal_flux(const std::vector<real_t> &state,
                                   const geo::Vec3 &n,
                                   std::vector<real_t> &normal_flux) const override;

    private:
        std::array<real_t, 4> compute_flux_x(real_t rho, real_t u, real_t v, real_t E, real_t p) const;
        std::array<real_t, 4> compute_flux_y(real_t rho, real_t u, real_t v, real_t E, real_t p) const;

        /// @brief Specific heat ratio
        real_t gamma_;
    };
}