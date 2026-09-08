#pragma once

#include <sfem/discretization/fvm/core/fv_field.hpp>
#include <filesystem>

namespace sfem::fvm
{    
    class FVSolver
    {
    public:
        FVSolver() = default;

        virtual void pre_timestep();

        virtual void post_timestep();

        virtual void step(real_t time, real_t dt) = 0;

    protected:
        /// @brief Current simulation time
        real_t time_ = 0.0;

        /// @brief Current timestep
        int step_ = 0;

        /// @brief Current timestep size
        real_t dt_ = 0.0;

        /// @brief Fields to plot
        std::vector<FVField> plot_fields_;

        /// @brief Root plot file
        std::filesystem::path plot_file_ = "solution";

        /// @brief Plot interval
        int plot_int_ = 10;
    };
}