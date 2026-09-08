#pragma once

#include <sfem/sfem.hpp>

using namespace sfem;
using namespace sfem::fvm;
using namespace sfem::fvm::simple;

SimpleOptions parse_simple_options(const nlohmann::json &j)
{
    SimpleOptions options;

    // Momentum equation options
    {
        SFEM_CHECK_JSON_ENTRY("momentum", j);
        const auto &momentum = j["momentum"];

        SFEM_CHECK_JSON_ENTRY("relaxation_factor", momentum);
        options.momentum_alpha = momentum["relaxation_factor"];

        if (momentum.contains("linear_solver"))
        {
            const auto &linear_solver = momentum["linear_solver"];

            SFEM_CHECK_JSON_ENTRY("type", linear_solver);
            options.momentum_solver_type = la::parse_linear_solver_type(linear_solver["type"]);

            SFEM_CHECK_JSON_ENTRY("options", linear_solver);
            options.momentum_solver_options = la::parse_linear_solver_options(linear_solver["options"]);
        }
    }

    // Pressure equation options
    {
        SFEM_CHECK_JSON_ENTRY("pressure", j);
        const auto &pressure = j["pressure"];

        SFEM_CHECK_JSON_ENTRY("relaxation_factor", pressure);
        options.momentum_alpha = pressure["relaxation_factor"];

        if (pressure.contains("linear_solver"))
        {
            const auto &linear_solver = pressure["linear_solver"];

            SFEM_CHECK_JSON_ENTRY("type", linear_solver);
            options.pressure_solver_type = la::parse_linear_solver_type(linear_solver["type"]);

            SFEM_CHECK_JSON_ENTRY("options", linear_solver);
            options.pressure_solver_options = la::parse_linear_solver_options(linear_solver["options"]);
        }
    }

    options.backend = la::parse_backend(j.value("backend", "native"));
    options.n_orthogonal_correctors = j.value("n_orthogonal_correctors", options.n_orthogonal_correctors);
    options.transient = j.value("transient", options.transient);
    options.atol_simple = j.value("atol", options.atol_simple);
    options.rtol_simple = j.value("rtol", options.rtol_simple);
    options.max_iter_simple = j.value("max_iter", options.max_iter_simple);
    options.plot_file = j.value("plot_file", options.plot_file);
    options.plot_int = j.value("plot_int", options.plot_int);

    return options;
}

SimpleSolver from_file(const std::filesystem::path &filename)
{
    // Open and parse JSON input file
    std::ifstream file(filename);
    SFEM_CHECK_FILE_OPEN(file, filename);
    const auto j = nlohmann::json::parse(file);

    // Read mesh
    SFEM_CHECK_JSON_ENTRY("mesh", j);
    auto mesh = io::read_mesh(j["mesh"], mesh::GhostMode::shared_facet);

    // Finite volume space
    auto V = std::make_shared<FVSpace>(mesh);

    SFEM_CHECK_JSON_ENTRY("fv_fields", j);
    const auto fv_fields = j["fv_fields"];

    // Velocity field
    std::vector<FVField> U;
    const std::array<std::string, 3> U_names = {"u", "v", "w"};
    for (int i = 0; i < mesh->pdim(); i++)
    {
        SFEM_CHECK_JSON_ENTRY(U_names[i], fv_fields);
        U.push_back(parse_field(fv_fields[U_names[i]], V, U_names[i]));
    }

    // Pressure field
    SFEM_CHECK_JSON_ENTRY("p", fv_fields);
    FVField p = parse_field( fv_fields["p"], V, "p");

    // Physical properties
    SFEM_CHECK_JSON_ENTRY("physical_properties", j);
    const auto &physical_properties = j["physical_properties"];

    // Density
    SFEM_CHECK_JSON_ENTRY("density", physical_properties);
    const real_t rho = physical_properties["density"];

    // Dynamic viscosity
    SFEM_CHECK_JSON_ENTRY("dynamic_viscosity", physical_properties);
    const real_t mu = physical_properties["dynamic_viscosity"];

    // SIMPLE options
    SFEM_CHECK_JSON_ENTRY("simple", j);
    SimpleOptions options = parse_simple_options(j["simple"]);

    return SimpleSolver(U, p, rho, mu, options);
}