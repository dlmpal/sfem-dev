#pragma once

#include "linear_solver.hpp"

namespace sfem::la
{
    enum class SolverType
    {
        gmres,
        cg
    };

    LinearSolver *create_solver(SolverType type, SolverOptions options);
}