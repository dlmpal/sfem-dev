#pragma once

#include <sfem/discretization/fvm/core/fv_solver.hpp>
#include <sfem/discretization/fvm/core/fv_equation.hpp>
#include <sfem/la/backend.hpp>

namespace sfem::fvm::simple
{
    struct SimpleOptions
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

        /// @brief Whether to include the transient term in the momentum equation
        bool transient = false;

        /// @brief SIMPLE absolute tolerance
        real_t atol_simple = 1e-10;

        /// @brief SIMPLE relative tolerance
        real_t rtol_simple = 1e-6;

        /// @brief Maximum number of SIMPLE iterations
        int max_iter_simple = 1000;

        /// @brief Root plot file
        std::filesystem::path plot_file = "solution";

        /// @brief Plot interval
        int plot_int = 10;
    };

    class SimpleSolver : public FVSolver
    {
    public:
        SimpleSolver(const std::vector<FVField> &U, const FVField &P,
                     real_t rho, real_t mu, SimpleOptions options);

        void pre_timestep() override;

        void step(real_t time, real_t dt) override;

    protected:
        void solve_momentum();

        void solve_pressure();

        void correct_fields();

        /// @brief Velocity field (per direction)
        std::vector<FVField> U_;

        /// @brief Old velocity field (per direction)
        std::vector<FVField> U_old_;

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
        std::vector<real_t> mdot_;

        /// @brief Momentum equation (per direction)
        std::vector<Equation> momentum_;

        /// @brief Pressure equation
        Equation pressure_;

        /// @brief Solver options
        SimpleOptions options_;

        /// @brief Mass imbalance of the last SIMPLE iteration
        real_t mass_residual_;

        /// @brief Residual norm history for the momentum and pressure
        /// equations, followed by the mass imbalance
        std::vector<real_t> residual_history_;
    };
}