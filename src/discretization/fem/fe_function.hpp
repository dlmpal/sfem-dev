#pragma once

#include "../function.hpp"
#include "fe_space.hpp"

namespace sfem::fem
{
    class FEFunction : public Function
    {
    public:
        FEFunction(std::shared_ptr<const FESpace> fe_space,
                   const std::vector<std::string> &components);

        std::shared_ptr<const FESpace> space() const;

        la::DenseMatrix cell_values(int cell_idx) const;

    private:
        std::shared_ptr<const FESpace> fe_space_;
    };
}