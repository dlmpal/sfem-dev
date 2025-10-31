#include "euler.hpp"
#include <cmath>

namespace sfem::fvm
{
    //=============================================================================
    EulerFlux::EulerFlux(real_t gamma, int dim)
        : FluxFunction(dim + 2, dim),
          gamma_(gamma)
    {
    }
    //=============================================================================
    real_t EulerFlux::compute_flux(const std::vector<real_t> &state,
                                   std::vector<real_t> &flux, int dir) const
    {
        // Get state variables
        const real_t rho = state[0];
        geo::Vec3 v;
        for (int i = 0; i < dim_; i++)
        {
            v(i) = state[i + 1] / rho;
        }
        const real_t v_mag = v.mag();
        const real_t KE = 0.5 * rho * v_mag * v_mag;
        const real_t E = state[dim_ + 1];
        const real_t e = (E - KE) / rho;
        const real_t p = (gamma_ - 1) * rho * e;

        // Compute fluxes
        flux[0] = rho * v(dir);
        for (int i = 0; i < dim_; i++)
        {
            flux[i + 1] = rho * v(i) * v(dir);
        }
        flux[dir + 1] += p;
        flux[dim_ + 1] = (E + p) * v(dir);

        // Max. wavespeed
        const double speed = v_mag;
        const double sound = std::sqrt(gamma_ * p / rho);

        return sound + speed;
    }
    //=============================================================================
    real_t EulerFlux::compute_normal_flux(const std::vector<real_t> &state,
                                          const geo::Vec3 &normal,
                                          std::vector<real_t> &normal_flux) const
    {
        // Get state variables
        const real_t rho = state[0];
        geo::Vec3 v;
        for (int i = 0; i < dim_; i++)
        {
            v(i) = state[i + 1] / rho;
        }
        const real_t v_mag = v.mag();
        const real_t KE = 0.5 * rho * v_mag * v_mag;
        const real_t E = state[dim_ + 1];
        const real_t e = (E - KE) / rho;
        const real_t p = (gamma_ - 1) * rho * e;

        // Compute normal speed
        const real_t v_normal = geo::inner(v, normal);

        // Compute fluxes
        normal_flux[0] = rho * v_normal;
        for (int i = 0; i < dim_; i++)
        {
            normal_flux[i + 1] = rho * v(i) * v_normal + p * normal(i);
        }
        normal_flux[dim_ + 1] = (E + p) * v_normal;

        // Max. wavespeed
        const double speed = v_mag;
        const double sound = std::sqrt(gamma_ * p / rho);

        return sound + speed;
    }
}