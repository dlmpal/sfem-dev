#pragma once

#include "linear_solver.hpp"

namespace sfem::la
{
    enum class SolverType
    {
        cg,
        gmres
    };

    LinearSolver *create_solver(SolverType type, SolverOptions options);
}