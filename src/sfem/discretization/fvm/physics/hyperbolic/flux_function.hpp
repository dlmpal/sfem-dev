#pragma once

#include <sfem/geo/vec3.hpp>
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
}