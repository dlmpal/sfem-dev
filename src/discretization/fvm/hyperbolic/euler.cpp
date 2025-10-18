#include "euler.hpp"
#include <cmath>

namespace sfem::fvm
{
    //=============================================================================
    Euler2D::Euler2D(real_t gamma)
        : FluxFunction(4, 2), gamma_(gamma)
    {
    }
    //=============================================================================
    std::array<real_t, 4> Euler2D::compute_flux_x(real_t rho, real_t u, real_t v, real_t E, real_t p) const
    {
        return {rho * u,
                rho * u * u + p,
                rho * u * v,
                u * (E + p)};
    }
    //=============================================================================
    std::array<real_t, 4> Euler2D::compute_flux_y(real_t rho, real_t u, real_t v, real_t E, real_t p) const
    {
        return {rho * v,
                rho * u * v,
                rho * v * v + p,
                v * (E + p)};
    }
    //=============================================================================
    real_t Euler2D::compute_flux(const std::vector<real_t> &state, std::vector<real_t> &flux, int dir) const
    {
        const real_t rho = state[0];
        const real_t u = state[1] / rho;
        const real_t v = state[2] / rho;
        const real_t E = state[3];
        const real_t KE = 0.5 * rho * (u * u + v * v);
        const real_t p = (gamma_ - 1.0) * (E - KE);

        std::array<real_t, 4> f;
        if (dir == 0)
        {
            f = compute_flux_x(rho, u, v, E, p);
        }
        else if (dir == 1)
        {
            f = compute_flux_y(rho, u, v, E, p);
        }

        flux[0] = f[0];
        flux[1] = f[1];
        flux[2] = f[2];
        flux[3] = f[3];

        const real_t sound = std::sqrt(gamma_ * p / rho);
        const real_t speed = std::sqrt(2.0 * KE / rho);

        return sound + speed;
    }
    //=============================================================================
    real_t Euler2D::compute_normal_flux(const std::vector<real_t> &state,
                                        const geo::Vec3 &n,
                                        std::vector<real_t> &normal_flux) const
    {
        const real_t rho = state[0];
        const real_t u = state[1] / rho;
        const real_t v = state[2] / rho;
        const real_t E = state[3];
        const real_t KE = 0.5 * rho * (u * u + v * v);
        const real_t p = (gamma_ - 1.0) * (E - KE);

        const auto f_x = compute_flux_x(rho, u, v, E, p);
        const auto f_y = compute_flux_y(rho, u, v, E, p);

        normal_flux[0] = f_x[0] * n.x() + f_y[0] * n.y();
        normal_flux[1] = f_x[1] * n.x() + f_y[1] * n.y();
        normal_flux[2] = f_x[2] * n.x() + f_y[2] * n.y();
        normal_flux[3] = f_x[3] * n.x() + f_y[3] * n.y();

        const real_t sound = std::sqrt(gamma_ * p / rho);
        const real_t speed = std::sqrt(2 * KE / rho);

        return sound + speed;
    }
}