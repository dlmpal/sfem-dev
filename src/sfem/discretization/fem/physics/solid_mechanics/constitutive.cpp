#include "constitutive.hpp"

namespace sfem::fem::solid_mechanics
{
    //=============================================================================
    LinearElasticIsotropic::LinearElasticIsotropic(Field &E, Field &nu, Field &rho)
        : E_(E),
          nu_(nu),
          rho_(rho)
    {
    }
    //=============================================================================
    Field &LinearElasticIsotropic::E()
    {
        return E_;
    }
    //=============================================================================
    const Field &LinearElasticIsotropic::E() const
    {
        return E_;
    }
    //=============================================================================
    Field &LinearElasticIsotropic::nu()
    {
        return nu_;
    }
    //=============================================================================
    const Field &LinearElasticIsotropic::nu() const
    {
        return nu_;
    }
    //=============================================================================
    Field &LinearElasticIsotropic::rho()
    {
        return rho_;
    }
    //=============================================================================
    const Field &LinearElasticIsotropic::rho() const
    {
        return rho_;
    }
    //=============================================================================
    LinearElasticPlaneStress::LinearElasticPlaneStress(Field &E, Field &nu, Field &rho)
        : LinearElasticIsotropic(E, nu, rho)
    {
    }
    //=============================================================================
    int LinearElasticPlaneStress::dim() const
    {
        return 2;
    }
    //=============================================================================
    int LinearElasticPlaneStress::n_strain() const
    {
        return 3;
    }
    //=============================================================================
    void LinearElasticPlaneStress::tangent(const FEData &data,
                                           const la::DenseMatrix &,
                                           la::DenseMatrix &D) const
    {
        const real_t E = E_.cell_value(data.elem_idx, data.pt);
        const real_t nu = nu_.cell_value(data.elem_idx, data.pt);
        const real_t c = E / (1 - nu * nu);
        D(0, 0) = c * 1.0;
        D(0, 1) = c * nu;
        D(1, 0) = c * nu;
        D(1, 1) = c * 1.0;
        D(2, 2) = c * (1 - nu) * 0.5;
    }
    //=============================================================================
    void LinearElasticPlaneStress::stress(const FEData &data,
                                          const la::DenseMatrix &e,
                                          la::DenseMatrix &s) const
    {
        la::DenseMatrix D(3, 3);
        tangent(data, e, D);
        s = D * e;
    }
    //=============================================================================
    LinearElastic3D::LinearElastic3D(Field &E, Field &nu, Field &rho)
        : LinearElasticIsotropic(E, nu, rho)
    {
    }
    //=============================================================================
    int LinearElastic3D::dim() const
    {
        return 3;
    }
    //=============================================================================
    int LinearElastic3D::n_strain() const
    {
        return 6;
    }
    //=============================================================================
    void LinearElastic3D::tangent(const FEData &data,
                                  const la::DenseMatrix &,
                                  la::DenseMatrix &D) const
    {
        const real_t E = E_.cell_value(data.elem_idx, data.pt);
        const real_t nu = nu_.cell_value(data.elem_idx, data.pt);
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
    }
    //=============================================================================
    void LinearElastic3D::stress(const FEData &data,
                                 const la::DenseMatrix &e,
                                 la::DenseMatrix &s) const
    {
        la::DenseMatrix D(6, 6);
        tangent(data, e, D);
        s = D * e;
    }
}