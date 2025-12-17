#include "linear_solver_factory.hpp"
#include "cg.hpp"
#include "gmres.hpp"

namespace sfem::la
{
    LinearSolver *create_solver(SolverType type, SolverOptions options)
    {
        LinearSolver *solver = nullptr;
        switch (type)
        {
        case SolverType::cg:
            solver = new CG(options);
            break;
        case SolverType::gmres:
            solver = new GMRES(options);
            break;
        default:
            break;
        }
        return solver;
    }
}