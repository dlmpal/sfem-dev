#include "numerical_flux.hpp"
#include <sfem/discretization/fvm/physics/hyperbolic/euler.hpp>

namespace sfem::fvm
{
    //=============================================================================
    NumericalFlux::NumericalFlux(std::shared_ptr<const FluxFunction> flux)
        : flux_(flux),
          flux1_(flux->n_comp()),
          flux2_(flux->n_comp())
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
    real_t RusanovFlux::compute_normal_flux(const std::vector<real_t> &state1,
                                            const std::vector<real_t> &state2,
                                            const geo::Vec3 &normal,
                                            std::vector<real_t> &normal_flux) const
    {
        // Compute left and right normal fluxes (and wave speeds)
        const real_t s1 = flux_->compute_normal_flux(state1, normal, flux1_);
        const real_t s2 = flux_->compute_normal_flux(state2, normal, flux2_);

        // Maximum wave speed
        const real_t s = std::max(std::abs(s1), std::abs(s2));

        // F = 0.5 * (F_R + F_L - |s| * (U_R - U_L))
        for (int i = 0; i < flux_->n_comp(); i++)
        {
            normal_flux[i] = 0.5 * (flux2_[i] + flux1_[i] - s * (state2[i] - state1[i]));
        }

        return s;
    }
    //=============================================================================
    GodunovFlux::GodunovFlux(std::shared_ptr<const FluxFunction> flux)
        : NumericalFlux(flux)
    {
    }
    //=============================================================================
    real_t GodunovFlux::compute_normal_flux(const std::vector<real_t> &state1,
                                            const std::vector<real_t> &state2,
                                            const geo::Vec3 &normal,
                                            std::vector<real_t> &normal_flux) const
    {
        // Compute left and right normal fluxes (and wave speeds)
        const real_t s1 = flux_->compute_normal_flux(state1, normal, flux1_);
        const real_t s2 = flux_->compute_normal_flux(state2, normal, flux2_);

        // Maximum wave speed
        const real_t s = std::max(std::abs(s1), std::abs(s2));

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            if (state1[i] <= state2[i])
            {
                normal_flux[i] = std::min(flux1_[i], flux2_[i]);
            }
            else
            {
                normal_flux[i] = std::max(flux1_[i], flux2_[i]);
            }
        }

        return s;
    }
    //=============================================================================
    HLLFlux::HLLFlux(std::shared_ptr<const EulerFlux> flux)
        : NumericalFlux(flux)
    {
    }
    //=============================================================================
    double HLLFlux::compute_normal_flux(const std::vector<real_t> &state1,
                                        const std::vector<real_t> &state2,
                                        const geo::Vec3 &normal,
                                        std::vector<real_t> &normal_flux) const
    {
        // Compute left and right normal fluxes (and wave speeds)
        const real_t s1_ = flux_->compute_normal_flux(state1, normal, flux1_);
        const real_t s2_ = flux_->compute_normal_flux(state2, normal, flux2_);

        // HLL wavespeeds
        const real_t s1 = -std::max(s1_, s2_);
        const real_t s2 = std::max(s1_, s2_);

        for (int i = 0; i < flux_->n_comp(); i++)
        {
            if (s1 > 0)
            {
                normal_flux[i] = flux1_[i];
            }
            else if (s2 > 0)
            {
                normal_flux[i] = (s2 * flux1_[i] - s1 * flux2_[i] + s1 * s2 * (state2[i] - state1[i])) / (s2 - s1);
            }
            else
            {
                normal_flux[i] = flux2_[i];
            }
        }

        return std::max(s1, s2);
    }
}