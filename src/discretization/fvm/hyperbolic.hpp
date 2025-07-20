#pragma once

#include "../../geo/vec3.hpp"
#include <vector>
#include <memory>

namespace sfem::fvm
{
    class FluxFunction
    {

    public:
        FluxFunction(int n_comp, int dim);

        int n_comp() const;

        int dim() const;

        /// @brief Compute the flux for a given direction
        virtual real_t compute_flux(const std::vector<real_t> &state,
                                    std::vector<real_t> &flux,
                                    int dir) const = 0;

        /// @brief Compute the inner product of the flux
        /// and a given normal vector.
        virtual real_t compute_normal_flux(const std::vector<real_t> &state,
                                           const geo::Vec3 &normal,
                                           std::vector<real_t> &normal_flux) const = 0;

    protected:
        /// @brief Number of components/equations
        int n_comp_;

        /// @brief Dimension
        int dim_;
    };

    class NumericalFlux
    {
    public:
        NumericalFlux(std::shared_ptr<const FluxFunction> flux);

        std::shared_ptr<const FluxFunction> flux_function() const;

        virtual real_t compute_normal_flux(const std::vector<real_t> &state_left,
                                           const std::vector<real_t> &state_right,
                                           const geo::Vec3 &normal,
                                           std::vector<real_t> &normal_flux) const = 0;

    protected:
        /// @brief Flux function
        std::shared_ptr<const FluxFunction> flux_;

        /// @brief Vectors to store left and right fluxes
        mutable std::vector<real_t> flux_left_, flux_right_;
    };

    class RusanovFlux : public NumericalFlux
    {
    public:
        RusanovFlux(std::shared_ptr<const FluxFunction> flux);

        real_t compute_normal_flux(const std::vector<real_t> &state_left,
                                   const std::vector<real_t> &state_right,
                                   const geo::Vec3 &normal,
                                   std::vector<real_t> &normal_flux) const override;
    };

    class GodunovFlux : public NumericalFlux
    {
    public:
        GodunovFlux(std::shared_ptr<const FluxFunction> flux);

        real_t compute_normal_flux(const std::vector<real_t> &state_left,
                                   const std::vector<real_t> &state_right,
                                   const geo::Vec3 &normal,
                                   std::vector<real_t> &normal_flux) const override;
    };
}