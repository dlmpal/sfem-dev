#include "linear_elasticity.hpp"

namespace sfem::fem::kernels
{
    //=============================================================================
    LinearElasticity2D::LinearElasticity2D(real_t E, real_t nu, real_t thick)
        : E_(E), nu_(nu), thick_(thick)
    {
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity2D::operator()(int cell_idx, const fem::FEData &data) const
    {
        // Stress-strain matrix
        la::DenseMatrix D(3, 3);
        real_t coeff = thick_ * E_ / (1 - nu_ * nu_);
        D(0, 0) = coeff * 1.0;
        D(0, 1) = coeff * nu_;
        D(1, 0) = coeff * nu_;
        D(1, 1) = coeff * 1.0;
        D(2, 2) = coeff * (1 - nu_) * 0.5;

        // Strain-displacement matrix
        const auto &dNdX = data.dNdX;
        int n_cols = dNdX.n_rows() * 2;
        int n_rows = 3;
        la::DenseMatrix B(n_rows, n_cols);
        for (int i = 0; i < dNdX.n_rows(); i++)
        {
            B(0, i * 2 + 0) = dNdX(i, 0); ///< exx
            B(1, i * 2 + 1) = dNdX(i, 1); ///< eyy
            B(2, i * 2 + 0) = dNdX(i, 1); ///< exy
            B(2, i * 2 + 1) = dNdX(i, 0); ///< exy
        }

        return B.T() * D * B;
    }
}