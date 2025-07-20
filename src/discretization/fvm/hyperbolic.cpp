#include "hyperbolic.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FluxFunction::FluxFunction(int n_comp, int dim)
        : n_comp_(n_comp),
          dim_(dim)
    {
    }
    //=============================================================================
    int FluxFunction::n_comp() const
    {
        return n_comp_;
    }
    //=============================================================================
    int FluxFunction::dim() const
    {
        return dim_;
    }
    //=============================================================================
    NumericalFlux::NumericalFlux(std::shared_ptr<const FluxFunction> flux)
        : flux_(flux),
          flux_left_(flux->n_comp()),
          flux_right_(flux->n_comp())
    {
    }
    //=============================================================================
    std::shared_ptr<const FluxFunction> NumericalFlux::flux_function() const
    {
        return flux_;
    }
    //=============================================================================
    RusanovFlux::RusanovFlux(std::shared_ptr<const FluxFunction> flux)
        : NumericalFlux(flux)
    {
    }
    //=============================================================================
    real_t RusanovFlux::compute_normal_flux(const std::vector<real_t> &state_left,
                                            const std::vector<real_t> &state_right,
                                            const geo::Vec3 &normal,
                                            std::vector<real_t> &normal_flux) const
    {
        // Compute left and right normal fluxes
        const real_t s_left = flux_->compute_normal_flux(state_left, normal, flux_left_);
        const real_t s_right = flux_->compute_normal_flux(state_right, normal, flux_right_);

        // Maximum wave speed
        const real_t s = std::max(std::abs(s_left), std::abs(s_right));

        // F = 0.5 * (F_R + F_L - |s| * (U_R - U_L))
        for (int i = 0; i < flux_->n_comp(); i++)
        {
            normal_flux[i] = 0.5 * (flux_right_[i] + flux_left_[i] - s * (state_right[i] - state_left[i]));
        }

        return s;
    }
    //=============================================================================
    GodunovFlux::GodunovFlux(std::shared_ptr<const FluxFunction> flux)
        : NumericalFlux(flux)
    {
    }
    //=============================================================================
    real_t GodunovFlux::compute_normal_flux(const std::vector<real_t> &state_left,
                                            const std::vector<real_t> &state_right,
                                            const geo::Vec3 &normal,
                                            std::vector<real_t> &normal_flux) const
    {
        // Compute left and right normal fluxes
        const real_t s_left = flux_->compute_normal_flux(state_left, normal, flux_left_);
        const real_t s_right = flux_->compute_normal_flux(state_right, normal, flux_right_);

        // Maximum wave speed
        const real_t s = std::max(std::abs(s_left), std::abs(s_right));

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            if (state_left[i] <= state_right[i])
            {
                normal_flux[i] = std::min(flux_left_[i], flux_right_[i]);
            }
            else
            {
                normal_flux[i] = std::max(flux_left_[i], flux_right_[i]);
            }
        }

        return s;
    }
}