#pragma once

#include "fv_space.hpp"
#include "../function.hpp"

namespace sfem::fvm
{
    /// @brief Function defined on a finite volume space
    class FVFunction : public Function
    {
    public:
        /// @brief Create a FVFunction
        /// @param fv_space Finite volume space
        /// @param components Component names
        FVFunction(std::shared_ptr<const FVSpace> fv_space,
                   const std::vector<std::string> &components);

        /// @brief Get the finite volume space
        std::shared_ptr<const FVSpace> space() const;

    private:
        /// @brief Finite volume space
        std::shared_ptr<const FVSpace> fv_space_;
    };
}