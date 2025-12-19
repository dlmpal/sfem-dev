#pragma once

#include <sfem/discretization/fvm/physics/hyperbolic/flux_function.hpp>

namespace sfem::fvm
{
    /// @brief Numerical flux ABC
    class NumericalFlux
    {
    public:
        /// @brief Create a NumericalFlux for a given FluxFunction
        NumericalFlux(std::shared_ptr<const FluxFunction> flux);

        /// @brief Get the flux function
        std::shared_ptr<const FluxFunction> flux_function() const;

        /// @brief Compute the normal flux at the interface of two cells/states
        virtual real_t compute_normal_flux(const std::vector<real_t> &state1,
                                           const std::vector<real_t> &state2,
                                           const geo::Vec3 &normal,
                                           std::vector<real_t> &normal_flux) const = 0;

    protected:
        /// @brief Flux function
        std::shared_ptr<const FluxFunction> flux_;

        /// @brief Vectors to store adjacent cell fluxes
        mutable std::vector<real_t> flux1_, flux2_;
    };

    /// @brief Rusanov numerical flux
    class RusanovFlux : public NumericalFlux
    {
    public:
        RusanovFlux(std::shared_ptr<const FluxFunction> flux);
        real_t compute_normal_flux(const std::vector<real_t> &state1,
                                   const std::vector<real_t> &state2,
                                   const geo::Vec3 &normal,
                                   std::vector<real_t> &normal_flux) const override;
    };

    /// @brief Godunov numerical flux
    class GodunovFlux : public NumericalFlux
    {
    public:
        GodunovFlux(std::shared_ptr<const FluxFunction> flux);
        real_t compute_normal_flux(const std::vector<real_t> &state1,
                                   const std::vector<real_t> &state2,
                                   const geo::Vec3 &normal,
                                   std::vector<real_t> &normal_flux) const override;
    };

    // Forward Declaration
    class EulerFlux;

    // Harten, Lax and Van Leer (HLL) approximate Riemann problem solver for the Euler equations
    class HLLFlux : public NumericalFlux
    {
    public:
        HLLFlux(std::shared_ptr<const EulerFlux> flux);
        real_t compute_normal_flux(const std::vector<real_t> &state1,
                                   const std::vector<real_t> &state2,
                                   const geo::Vec3 &normal,
                                   std::vector<real_t> &normal_flux) const override;
    };
}