#include "linear_elasticity.hpp"

namespace sfem::fem::kernels
{
    //=============================================================================
    LinearElasticity2D::LinearElasticity2D(std::shared_ptr<Function> E,
                                           std::shared_ptr<Function> nu,
                                           std::shared_ptr<Function> thick)
        : E_(E), nu_(nu), thick_(thick)
    {
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity2D::operator()(int cell_idx, const fem::FEData &data) const
    {
        // Cell coefficients
        const real_t E = (*E_)(cell_idx, 0);
        const real_t nu = (*nu_)(cell_idx, 0);
        const real_t thick = (*thick_)(cell_idx, 0);

        // Stress-strain matrix
        la::DenseMatrix D(3, 3);
        const real_t coeff = thick * E / (1 - nu * nu);
        D(0, 0) = coeff * 1.0;
        D(0, 1) = coeff * nu;
        D(1, 0) = coeff * nu;
        D(1, 1) = coeff * 1.0;
        D(2, 2) = coeff * (1 - nu) * 0.5;

        // Strain-displacement matrix
        const auto &dNdX = data.dNdX;
        const int n_cols = dNdX.n_rows() * 2;
        const int n_rows = 3;
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
    //=============================================================================
    LinearElasticity3D::LinearElasticity3D(std::shared_ptr<Function> E,
                                           std::shared_ptr<Function> nu)
        : E_(E), nu_(nu)
    {
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity3D::operator()(int cell_idx, const fem::FEData &data) const
    {
        // Cell coefficients
        const real_t E = (*E_)(cell_idx, 0);
        const real_t nu = (*nu_)(cell_idx, 0);

        // Stress-strain matrix
        la::DenseMatrix D(6, 6);
        const real_t c1 = E / ((1 + nu) * (1 - 2 * nu));
        const real_t c2 = (1 - 2 * nu) / 2;
        D(0, 0) = (1 - nu) * c1;
        D(0, 1) = nu * c1;
        D(0, 2) = nu * c1;
        D(1, 0) = nu * c1;
        D(1, 1) = (1 - nu) * c1;
        D(1, 2) = nu * c1;
        D(2, 0) = nu * c1;
        D(2, 1) = nu * c1;
        D(2, 2) = (1 - nu) * c1;
        D(3, 3) = c1 * c2;
        D(4, 4) = c1 * c2;
        D(5, 5) = c1 * c2;

        // Strain-displacement matrix
        const auto &dNdX = data.dNdX;
        const int n_cols = dNdX.n_rows() * 3;
        const int n_rows = 6;
        la::DenseMatrix B(n_rows, n_cols);
        for (int i = 0; i < dNdX.n_rows(); i++)
        {
            B(0, i * 3 + 0) = dNdX(i, 0); ///< exx
            B(1, i * 3 + 1) = dNdX(i, 1); ///< eyy
            B(2, i * 3 + 2) = dNdX(i, 2); ///< ezz

            B(3, i * 3 + 0) = dNdX(i, 1); ///< exy
            B(3, i * 3 + 1) = dNdX(i, 0); ///< exy

            B(4, i * 3 + 1) = dNdX(i, 2); ///< eyz
            B(4, i * 3 + 2) = dNdX(i, 1); ///< eyz

            B(5, i * 3 + 0) = dNdX(i, 2); ///< ezx
            B(5, i * 3 + 2) = dNdX(i, 0); ///< ezx
        }

        return B.T() * D * B;
    }
}