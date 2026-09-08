#include "fv_solver.hpp"

#include <sfem/io/vtk/vtk.hpp>
#include <format>

namespace sfem::fvm
{
    //=============================================================================
    void FVSolver::pre_timestep()
    {
    }
    //=============================================================================
    void FVSolver::post_timestep()
    {
        if (!plot_fields_.empty() && step_ % plot_int_ == 0)
        {
            io::vtk::write(std::format("{}_{}", plot_file_.c_str(), step_), plot_fields_);
        }
    }
}