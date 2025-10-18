#pragma once

#include "fv_space.hpp"
#include "../field.hpp"

namespace sfem::fvm
{
    /// @brief Field defined on a finite volume space
    class FVField : public Field
    {
    public:
        /// @brief Create a FVFunction
        /// @param fv_space Finite volume space
        /// @param components Component names
        FVField(std::shared_ptr<const FVSpace> fv_space,
                const std::vector<std::string> &components);

        /// @brief Get the finite volume space
        std::shared_ptr<const FVSpace> space() const;

    private:
        /// @brief Finite volume space
        std::shared_ptr<const FVSpace> fv_space_;
    };
}