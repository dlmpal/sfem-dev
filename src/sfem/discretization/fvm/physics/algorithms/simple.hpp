#pragma once

#include <sfem/discretization/fvm/core/fv_equation.hpp>
#include <sfem/la/backend.hpp>

namespace sfem::fvm::algo
{
    struct SIMPLEOptions
    {
        /// @brief Momentum under-relaxation factor
        real_t momentum_alpha = 0.7;

        /// @brief Pressure under-relaxation factor
        real_t pressure_alpha = 0.3;

        /// @brief Momentum linear solver type
        la::SolverType momentum_solver_type = la::SolverType::gmres;

        /// @brief Momentum linear solver options
        la::SolverOptions momentum_solver_options = {};

        /// @brief Pressure linear solver type
        la::SolverType pressure_solver_type = la::SolverType::cg;

        /// @brief Pressure linear solver options
        la::SolverOptions pressure_solver_options = {};

        /// @brief Linear algebra backend
        la::Backend backend = la::Backend::native;

        /// @brief Number of orthogonal-correction iterations
        int n_orthogonal_correctors = 0;

        /// @brief Whether the flow is transient
        bool transient = false;

        /// @brief Maximum number of SIMPLE iterations
        int max_iter_simple = 50;

        /// @brief SIMPLE relative tolerance
        real_t rtol_simple = 1e-4;

        /// @brief Plot interval
        int plot_int = 10;
    };

    class SIMPLESolver
    {
    public:
        SIMPLESolver(std::vector<FVField> U, FVField P,
                     real_t rho, real_t mu, SIMPLEOptions options);

        void step(real_t time, real_t dt);

    private:
        void assemble_pressure_diffusivity();
        void assemble_mass_flux();
        void correct_fields();

    protected:
        /// @brief Velocity field (per direction)
        std::vector<FVField> U_;

        /// @brief Pressure field
        FVField P_;

        /// @brief Pressure correction field
        FVField Pcorr_;

        /// @brief Pressure diffusivity
        FVField D_;

        /// @brief Density
        ConstantField rho_;

        /// @brief Dynamic viscosity
        ConstantField mu_;

        /// @brief Mass flux
        /// @todo
        std::vector<real_t> flux_;

        /// @brief Momentum equation (per direction)
        std::vector<Equation> momentum_;

        /// @brief Pressure equation
        Equation pressure_;

        /// @brief Solver options
        SIMPLEOptions options_;

        /// @brief Current timestep
        real_t dt_;
    };
}