#include "linear_solver_factory.hpp"
#include <sfem/la/native/linear_solvers/gmres.hpp>
#include <sfem/la/native/linear_solvers/cg.hpp>

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