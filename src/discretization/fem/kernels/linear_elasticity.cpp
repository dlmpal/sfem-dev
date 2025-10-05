#include "linear_elasticity.hpp"
#include <cmath>

namespace sfem::fem::kernels
{
    //=============================================================================
    LinearElasticity2D::LinearElasticity2D(std::shared_ptr<const Coefficient> E,
                                           std::shared_ptr<const Coefficient> nu,
                                           std::shared_ptr<const Coefficient> thick)
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
    InertialLoad2D::InertialLoad2D(std::shared_ptr<const Coefficient> thick,
                                   std::shared_ptr<const Coefficient> rho,
                                   const std::array<real_t, 2> &g)
        : thick_(thick), rho_(rho), g_(g)
    {
    }
    //=============================================================================
    la::DenseMatrix InertialLoad2D::operator()(int cell_idx, const fem::FEData &data)
    {
        const real_t thick = (*thick_)(cell_idx, 0);
        const real_t rho = (*rho_)(cell_idx, 0);
        la::DenseMatrix F(data.N.n_rows() * 2, 1);
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            F(i * 2 + 0, 0) = rho * thick * g_[0] * data.N(i, 0);
            F(i * 2 + 1, 0) = rho * thick * g_[1] * data.N(i, 0);
        }
        return F;
    }
    //=============================================================================
    PressureLoad2D::PressureLoad2D(std::shared_ptr<const Coefficient> thick,
                                   std::shared_ptr<const Coefficient> pressure)
        : thick_(thick), pressure_(pressure)
    {
    }
    //=============================================================================
    la::DenseMatrix PressureLoad2D::operator()(int facet_idx,
                                               const fem::FEData &data,
                                               const geo::Vec3 &normal)
    {
        const real_t thick = (*thick_)(facet_idx, 0);
        const real_t pressure_value = (*pressure_)(facet_idx, 0);
        la::DenseMatrix F(data.N.n_rows() * 2, 1, 0.0);
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            F(i * 2 + 0, 0) += -pressure_value * thick * data.N(i, 0) * normal.x();
            F(i * 2 + 1, 0) += -pressure_value * thick * data.N(i, 0) * normal.y();
        }
        return F;
    }
    //=============================================================================
    la::DenseMatrix stress_strain3D(real_t E, real_t nu)
    {
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
        return D;
    }
    //=============================================================================
    la::DenseMatrix strain_displacement3D(const la::DenseMatrix &dNdX)
    {
        // Strain-displacement matrix
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
        return B;
    }
    //=============================================================================
    LinearElasticity3D::LinearElasticity3D(std::shared_ptr<const Coefficient> E,
                                           std::shared_ptr<const Coefficient> nu)
        : E_(E), nu_(nu)
    {
    }
    //=============================================================================
    la::DenseMatrix LinearElasticity3D::operator()(int cell_idx, const fem::FEData &data) const
    {
        // Cell coefficients
        const real_t E = (*E_)(cell_idx, 0);
        const real_t nu = (*nu_)(cell_idx, 0);

        // Stress-strain and strain-displacement matrices
        const auto D = stress_strain3D(E, nu);
        const auto B = strain_displacement3D(data.dNdX);

        // Stiffness matrix
        return B.T() * D * B;
    }
    //=============================================================================
    VonMises3D::VonMises3D(std::shared_ptr<const Coefficient> E,
                           std::shared_ptr<const Coefficient> nu,
                           std::shared_ptr<const FEFunction> U)
        : E_(E), nu_(nu), U_(U)
    {
    }
    //=============================================================================
    la::DenseMatrix VonMises3D::operator()(int cell_idx, const fem::FEData &data) const
    {
        // Cell coefficients
        const real_t E = (*E_)(cell_idx, 0);
        const real_t nu = (*nu_)(cell_idx, 0);

        // Stress-strain and strain-displacement matrices
        const auto D = stress_strain3D(E, nu);
        const auto B = strain_displacement3D(data.dNdX);

        // Nodal displacements
        const auto U = U_->cell_values(cell_idx);

        // Stresses
        const auto S = D * B * U;
        const real_t s_xx = S(0, 0);
        const real_t s_yy = S(1, 0);
        const real_t s_zz = S(2, 0);
        const real_t s_xy = S(3, 0);
        const real_t s_yz = S(4, 0);
        const real_t s_xz = S(5, 0);

        // Von-Mises stress
        real_t s_vm = std::pow(s_xx - s_yy, 2) + std::pow(s_yy - s_zz, 2) + std::pow(s_zz - s_xx, 2);
        s_vm += 6 * (std::pow(s_xy, 2) + std::pow(s_yz, 2) + std::pow(s_xz, 2));
        s_vm = std::sqrt(0.5 * s_vm);

        return la::DenseMatrix(1, 1, s_vm);
    }
    //=============================================================================
    InertialLoad3D::InertialLoad3D(std::shared_ptr<const Coefficient> rho,
                                   const std::array<real_t, 3> &g)
        : rho_(rho), g_(g)
    {
    }
    //=============================================================================
    la::DenseMatrix InertialLoad3D::operator()(int cell_idx, const fem::FEData &data)
    {
        const real_t rho = (*rho_)(cell_idx, 0);
        la::DenseMatrix F(data.N.n_rows() * 3, 1);
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            F(i * 3 + 0, 0) = rho * g_[0] * data.N(i, 0);
            F(i * 3 + 1, 0) = rho * g_[1] * data.N(i, 0);
            F(i * 3 + 2, 0) = rho * g_[2] * data.N(i, 0);
        }
        return F;
    }
    //=============================================================================
    PressureLoad3D::PressureLoad3D(std::shared_ptr<const Coefficient> pressure)
        : pressure_(pressure)
    {
    }
    //=============================================================================
    la::DenseMatrix PressureLoad3D::operator()(int facet_idx,
                                               const fem::FEData &data,
                                               const geo::Vec3 &normal)
    {
        const real_t pressure_value = (*pressure_)(facet_idx, 0);
        la::DenseMatrix F(data.N.n_rows() * 3, 1);
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            F(i * 2 + 0, 0) += -pressure_value * data.N(i, 0) * normal.x();
            F(i * 2 + 1, 0) += -pressure_value * data.N(i, 0) * normal.y();
            F(i * 3 + 2, 0) += -pressure_value * data.N(i, 0) * normal.z();
        }
        return F;
    }
}