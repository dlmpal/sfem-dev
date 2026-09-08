#include "la.hpp"
#include <sfem/base/error.hpp>

namespace sfem::la
{
    //=============================================================================
    la::Backend parse_backend(const std::string &backend_str)
    {
        la::Backend backend;
        if (backend_str == "native")
        {
            backend = la::Backend::native;
        }
        else if (backend_str == "petsc")
        {
            backend = la::Backend::petsc;
        }
        else
        {
            SFEM_ERROR(std::format("Unsupported linear algebra backend: {}\n", backend_str));
        }
        return backend;
    }
    //=============================================================================
    la::SolverType parse_linear_solver_type(const std::string &type_str)
    {
        la::SolverType type;
        if (type_str == "gmres")
        {
            type = la::SolverType::gmres;
        }
        else if (type_str == "cg")
        {
            type = la::SolverType::cg;
        }
        else
        {
            SFEM_ERROR(std::format("Unsupported linear solver type: {}\n", type_str));
        }
        return type;
    }
    //=============================================================================
    la::SolverOptions parse_linear_solver_options(const nlohmann::json &j)
    {
        la::SolverOptions options;
        options.atol = j.value("atol", options.atol);
        options.rtol = j.value("rtol", options.rtol);
        options.dtol = j.value("dtol", options.dtol);
        options.n_iter_max = j.value("n_iter_max", options.n_iter_max);
        options.print_conv = j.value("print_conv", options.print_conv);
        options.print_iter = j.value("print_iter", options.print_iter);
        return options;
    }
}