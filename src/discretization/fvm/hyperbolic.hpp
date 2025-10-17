#pragma once

#include "../../geo/vec3.hpp"
#include <vector>
#include <memory>

namespace sfem::fvm
{
    /// @brief Flux function ABC
    class FluxFunction
    {

    public:
        /// @brief Create a flux function
        /// @param n_comp Number of components/equations
        /// @param dim Dimension
        FluxFunction(int n_comp, int dim);

        /// @brief Get the number of components (i.e. equations)
        int n_comp() const;

        /// @brief Get the flux's dimensionality
        int dim() const;

        /// @brief Compute the flux for a given direction
        virtual real_t compute_flux(const std::vector<real_t> &state,
                                    std::vector<real_t> &flux, int dir) const = 0;

        /// @brief Compute the inner product of the flux and a given normal vector.
        virtual real_t compute_normal_flux(const std::vector<real_t> &state,
                                           const geo::Vec3 &normal,
                                           std::vector<real_t> &normal_flux) const = 0;

    protected:
        /// @brief Number of components/equations
        int n_comp_;

        /// @brief Dimension
        int dim_;
    };

    /// @brief Flux function for the two-dimensional Euler equations
    /// @note Subject to change
    class Euler2D : public FluxFunction
    {
    public:
        Euler2D(real_t gamma);

        real_t compute_flux(const std::vector<real_t> &state, std::vector<real_t> &flux, int dir) const;

        real_t compute_normal_flux(const std::vector<real_t> &state,
                                   const geo::Vec3 &n,
                                   std::vector<real_t> &normal_flux) const;

    private:
        std::array<real_t, 4> compute_flux_x(real_t rho, real_t u, real_t v, real_t E, real_t p) const;
        std::array<real_t, 4> compute_flux_y(real_t rho, real_t u, real_t v, real_t E, real_t p) const;

        /// @brief Specific heat ratio
        real_t gamma_;
    };

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
}