#pragma once

#include "quadrature.hpp"

namespace sfem::fem::quadrature
{
    template <std::size_t dim>
    class Gauss : public IntegrationRule
    {
    public:
        Gauss(int order);

        /// @brief Set the number of integration points per direction
        void set_n_points(int n_points) override;

        IntegrationPoint point(int i) const override;

    private:
        int order_;
    };
}