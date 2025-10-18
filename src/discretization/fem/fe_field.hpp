#pragma once

#include "../field.hpp"
#include "fe_space.hpp"

namespace sfem::fem
{
    /// @brief Field defined on a finite element space
    class FEField : public Field
    {
    public:
        /// @brief Create a FEField for a given space
        /// @param fe_space Finite element space
        /// @param components Component names
        FEField(std::shared_ptr<const FESpace> fe_space,
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