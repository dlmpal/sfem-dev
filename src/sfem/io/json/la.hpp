#pragma once

#include <sfem/io/json/json.hpp>
#include <sfem/la/backend.hpp>
#include <sfem/la/native/linear_solvers/linear_solver_factory.hpp>

namespace sfem::la
{
    la::Backend parse_backend(const std::string &backend_str);
    la::SolverType parse_linear_solver_type(const std::string &type_str);
    la::SolverOptions parse_linear_solver_options(const nlohmann::json &j);
}