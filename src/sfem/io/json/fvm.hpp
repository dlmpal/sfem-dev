#pragma once

#include <sfem/io/json/json.hpp>
#include <sfem/discretization/fvm/core/fv_field.hpp>

namespace sfem::fvm
{
    FVField parse_field(const nlohmann::json &j,
                        std::shared_ptr<const FVSpace> V,
                        const std::string &name);
}