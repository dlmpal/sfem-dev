#include "fvm.hpp"
#include <sfem/la/native/vector.hpp>

namespace sfem::fvm
{
    //=============================================================================
    FVField parse_field(const nlohmann::json &j,
                        std::shared_ptr<const FVSpace> V,
                        const std::string &name)
    {
        // Get gradient method
        SFEM_CHECK_JSON_ENTRY("gradient_method", j);
        GradientMethod gradient_method = GradientMethod::none;
        if (j["gradient_method"] == "green_gauss")
        {
            gradient_method = GradientMethod::green_gauss;
        }
        else if (j["gradient_method"] == "least_squares")
        {
            gradient_method = GradientMethod::least_squares;
        }

        // Create the field
        FVField phi(V, {name}, gradient_method);

        // Set boundary conditions
        SFEM_CHECK_JSON_ENTRY("boundary_conditions", j);
        {
            FVBC &phi_bc = phi.boundary_condition();
            for (const auto &[region_name, bc] : j["boundary_conditions"].items())
            {
                SFEM_CHECK_JSON_ENTRY("type", bc);
                BCType type = BCType::zero_neumann;
                if (bc["type"] == "neumann")
                {
                    type = BCType::neumann;
                }
                else if (bc["type"] == "dirichlet")
                {
                    type = BCType::dirichlet;
                }
                else if (bc["type"] == "robin")
                {
                    type = BCType::robin;
                }

                BCData data;
                if (type != BCType::zero_neumann)
                {
                    if (bc.contains("value"))
                    {
                        data.c = bc["value"];
                    }
                    else
                    {
                        SFEM_CHECK_JSON_ENTRY("a", bc);
                        SFEM_CHECK_JSON_ENTRY("b", bc);
                        SFEM_CHECK_JSON_ENTRY("c", bc);
                        data.a = bc["a"];
                        data.b = bc["b"];
                        data.c = bc["c"];
                    }
                }

                phi_bc.set_region_bc(region_name, type, data);
            }
        }

        // Set initial condition
        SFEM_CHECK_JSON_ENTRY("initial_condition", j);
        {
            /// @todo Implement ICType
            // ICType type = ICType::uniform;
            // if (type == ICType::uniform)
            if (true)
            {
                SFEM_CHECK_JSON_ENTRY("value", j["initial_condition"]);
                const real_t value = j["initial_condition"]["value"];
                phi.values().set_all(value);
            }
            else
            {
                SFEM_ERROR("Only uniform FV field initilization supported currently\n");
            }

            // Update the gradient after initialization
            phi.values().update_ghosts();
            phi.update_gradient();
        }

        return phi;
    }
}