#pragma once

#include "../function.hpp"
#include "fe_space.hpp"

namespace sfem::fem
{
    /// @brief Function defined on a finite element space
    class FEFunction : public Function
    {
    public:
        /// @brief Create a FEFunction for a given space
        /// @param fe_space Finite element space
        /// @param components Component names
        FEFunction(std::shared_ptr<const FESpace> fe_space,
                   const std::vector<std::string> &components);

        /// @brief Get the finite element space
        std::shared_ptr<const FESpace> space() const;

        /// @brief Get the values for a given cell
        la::DenseMatrix cell_values(int cell_idx) const;

    private:
        /// @brief Finite element space
        std::shared_ptr<const FESpace> fe_space_;
    };
}