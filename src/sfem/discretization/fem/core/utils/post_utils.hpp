#pragma once

#include <sfem/discretization/fem/core/fe_field.hpp>

namespace sfem::fem
{
    /// @brief Evaluate the quadrature point average value for each cell
    /// for a given operator
    /// @note F should be a cell constant field, i.e. have order = 0
    void cell_qpoint_average(const FESpace &V, ElementOperator op, FEField &F);
}