#pragma once

#include "fv_space.hpp"
#include "../function.hpp"

namespace sfem::fvm
{
    class FVFunction : public Function
    {
    public:
        FVFunction(std::shared_ptr<const FVSpace> fv_space,
                   const std::vector<std::string> &components);

        std::shared_ptr<const FVSpace> space() const;

    private:
        std::shared_ptr<const FVSpace> fv_space_;
    };
}