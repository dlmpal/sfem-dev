#pragma once

#include "fv_function.hpp"

namespace sfem::fvm
{
    enum class BCType
    {
        zero_neumann,
        neumann,
        dirichlet
    };

    class FVBC
    {
    public:
        FVBC(std::shared_ptr<const FVFunction> phi);

        // Avoid copies
        FVBC(const FVBC &) = delete;
        FVBC &operator=(const FVBC &) = delete;

        // Move constructor and assignment
        FVBC(FVBC &&) = default;
        FVBC &operator=(FVBC &&) = default;

        BCType region_type(const std::string &region) const;

        void set_value(const std::string &region,
                       const std::string &comp_name,
                       BCType type, real_t value);

        std::shared_ptr<const FVFunction> phi_;
        std::unordered_map<std::string, BCType> types_;
        std::unordered_map<int, real_t> values_;
    };
}