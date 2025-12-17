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

        real_t weight(int qpt_idx) const override;
        std::array<real_t, 3> point(int qpt_idx) const override;

    private:
        int order_;
    };
}