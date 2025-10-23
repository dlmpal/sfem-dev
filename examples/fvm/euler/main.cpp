// Solve the linearized 2D Euler equations for an acoustic monopole source
// Create the mesh with the cart-mesh application by:
// ${SFEM_DEV_INSTALL_DIR}/bin/cart-mesh -d=2 -Nx=200 -Ny=200 -x-low=-100 -x-high=100 -y-low=-100 -y-high=100

#include "sfem.hpp"

using namespace sfem;

namespace sfem::fvm
{
    struct AccousticParams
    {
        AccousticParams()
        {
            sound = std::sqrt(temp * gamma * R);
        }

        real_t rho = 1.293;
        real_t temp = 289;

        real_t gamma = 1.4;
        real_t R = 287.0;

        real_t freq = 14;
        real_t alpha = 0.167;
        real_t mach = 0.0;
        real_t sound;

        std::array<real_t, 3> center = {};
    };

    class Monopole
    {
    public:
        Monopole(const AccousticParams &params)
            : params_(params)
        {
        }

        void operator()(const std::array<real_t, 3> &pt,
                        std::span<real_t> values,
                        real_t time)
        {
            const real_t r = geo::compute_distance(pt, params_.center);
            const real_t c = params_.sound;
            const real_t f = params_.freq;
            const real_t a = params_.alpha;
            values[0] = std::sin(2 * M_PI * f * time) * std::exp(-a * r * r) / (c / c);
        }

    private:
        AccousticParams params_;
    };

    std::array<real_t, 3> compute_flux_x(real_t rho0, real_t c0, real_t u0,
                                         real_t rho, real_t u, real_t v)
    {
        return {u0 * rho + rho0 * u,
                c0 * c0 / rho0 * rho + u0 * u,
                u0 * v};
    }

    std::array<real_t, 3> compute_flux_y(real_t rho0, real_t c0, real_t u0,
                                         real_t rho, real_t u, real_t v)
    {
        return {rho0 * v,
                0.0,
                c0 * c0 / rho0 * rho};
    }

    class LinearizedEuler2D : public FluxFunction
    {

    public:
        LinearizedEuler2D(const AccousticParams &params)
            : FluxFunction(3, 2), params_(params)
        {
        }

        real_t compute_flux(const std::vector<real_t> &state,
                            std::vector<real_t> &flux,
                            int dir) const
        {
            const real_t rho0 = params_.rho;
            const real_t c0 = params_.sound;
            const real_t u0 = params_.mach * c0;

            const real_t rho = state[0];
            const real_t u = state[1];
            const real_t v = state[2];

            std::array<real_t, 3> f;
            if (dir == 0)
            {
                f = compute_flux_x(rho0, c0, u0,
                                   rho, u, v);
            }
            else if (dir == 1)
            {
                f = compute_flux_y(rho0, c0, u0,
                                   rho, u, v);
            }

            flux[0] = f[0];
            flux[1] = f[1];
            flux[2] = f[2];

            return c0 + u0;
        }

        real_t compute_normal_flux(const std::vector<real_t> &state,
                                   const geo::Vec3 &n,
                                   std::vector<real_t> &normal_flux) const
        {
            const real_t rho0 = params_.rho;
            const real_t c0 = params_.sound;
            const real_t u0 = params_.mach * c0;

            const real_t rho = state[0];
            const real_t u = state[1];
            const real_t v = state[2];

            const auto f_x = compute_flux_x(rho0, c0, u0,
                                            rho, u, v);
            const auto f_y = compute_flux_y(rho0, c0, u0,
                                            rho, u, v);

            normal_flux[0] = f_x[0] * n.x() + f_y[0] * n.y();
            normal_flux[1] = f_x[1] * n.x() + f_y[1] * n.y();
            normal_flux[2] = f_x[2] * n.x() + f_y[2] * n.y();

            return c0 + u0;
        }

    private:
        AccousticParams params_;
    };
}

int main(int argc, char **argv)
{
    initialize(argc, argv);

    auto mesh = io::read_mesh(argv[1], mesh::PartitionCriterion::shared_facet);

    auto V = std::make_shared<fvm::FVSpace>(mesh);
    auto S_new = std::make_shared<fvm::FVField>(V, std::vector<std::string>{"rho", "u", "v"});
    auto S_old = std::make_shared<fvm::FVField>(V, std::vector<std::string>{"rho", "u", "v"});
    auto rhs = std::make_shared<fvm::FVField>(V, std::vector<std::string>{"rho", "u", "v"});

    const real_t dt = 1e-3;
    const real_t time_start = 0.;
    const real_t time_stop = 0.1;
    const int plot_int = 5;

    int step = 0;
    real_t time = 0.0;

    fvm::AccousticParams params;
    fvm::Monopole src(params);

    // U_(n+1) = U_n + dt (Q - dV * sum_([F(UL,UR), G(UL,UR)] * [n_x, n_y]))
    auto flux = std::make_shared<fvm::LinearizedEuler2D>(params);
    auto nflux = std::make_shared<fvm::RusanovFlux>(flux);
    auto integrator = ode::create_erk(*S_new, fvm::ode::create_rhs(S_new, nflux, src), ode::ERKType::rk4);

    while (time <= time_stop)
    {
        // Increment time
        step++;
        time += dt;
        log_msg(std::format("Time: {}, Step: {}\n", time, step), true);

        // Advance
        integrator.advance(*S_old, *S_new, time, dt);

        // Save solution to file
        if (step % plot_int == 0)
        {
            S_new->update_ghosts();
            io::vtk::write(std::format("post/solution_{}", step), *mesh, {S_new});
        }

        la::copy(*S_new, *S_old);
    }

    return 0;
}