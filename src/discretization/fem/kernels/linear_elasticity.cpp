#include "linear_elasticity.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fem::kernels::elasticity
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
    LinearElasticPlaneStress::LinearElasticPlaneStress(Field &E, Field &nu,
                                                       Field &rho, Field &thick)
        : LinearElasticIsotropic(E, nu, rho),
          thick_(thick)
    {
    }
    //=============================================================================
    Field &LinearElasticPlaneStress::thick()
    {
        return thick_;
    }
    //=============================================================================
    const Field &LinearElasticPlaneStress::thick() const
    {
        return thick_;
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
    void LinearElasticPlaneStress::tangent(int cell_idx,
                                           mesh::CellType cell_type,
                                           const std::array<real_t, 3> &pt,
                                           la::DenseMatrix &D) const
    {
        const real_t E = E_.value(cell_idx, cell_type, pt);
        const real_t nu = nu_.value(cell_idx, cell_type, pt);
        const real_t t = thick_.value(cell_idx, cell_type, pt);
        const real_t c = t * E / (1 - nu * nu);
        D(0, 0) = c * 1.0;
        D(0, 1) = c * nu;
        D(1, 0) = c * nu;
        D(1, 1) = c * 1.0;
        D(2, 2) = c * (1 - nu) * 0.5;
    }
    //=============================================================================
    void LinearElasticPlaneStress::stress(int cell_idx,
                                          mesh::CellType cell_type,
                                          const std::array<real_t, 3> &pt,
                                          const la::DenseMatrix &strain,
                                          la::DenseMatrix &stress) const
    {
        la::DenseMatrix D(3, 3);
        tangent(cell_idx, cell_type, pt, D);
        stress = D * strain;
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
    void LinearElastic3D::tangent(int cell_idx,
                                  mesh::CellType cell_type,
                                  const std::array<real_t, 3> &pt,
                                  la::DenseMatrix &D) const
    {
        const real_t E = E_.value(cell_idx, cell_type, pt);
        const real_t nu = nu_.value(cell_idx, cell_type, pt);
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
    void LinearElastic3D::stress(int cell_idx,
                                 mesh::CellType cell_type,
                                 const std::array<real_t, 3> &pt,
                                 const la::DenseMatrix &strain,
                                 la::DenseMatrix &stress) const
    {
        la::DenseMatrix D(6, 6);
        tangent(cell_idx, cell_type, pt, D);
        stress = D * strain;
    }
    //=============================================================================
    LinearElasticity::LinearElasticity(FEField U,
                                       LinearElasticIsotropic &constitutive,
                                       const std::array<real_t, 3> &g)
        : U_(U),
          constitutive_(constitutive),
          g_(g)
    {
        if (U_.n_comp() != U_.space()->mesh()->pdim() ||
            U_.n_comp() != constitutive_.dim())
        {
            SFEM_ERROR("Mismatch between displacement field, constitutive law and mesh dimensions\n");
        }
    }
    //=============================================================================
    void LinearElasticity::strain_displacement(const FEData &data, la::DenseMatrix &B) const
    {
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            if (constitutive_.dim() == 2)
            {
                B(0, i * 2 + 0) = data.dNdX(i, 0); ///< exx
                B(1, i * 2 + 1) = data.dNdX(i, 1); ///< eyy
                B(2, i * 2 + 0) = data.dNdX(i, 1); ///< exy
                B(2, i * 2 + 1) = data.dNdX(i, 0); ///< exy
            }
            else
            {
                B(0, i * 3 + 0) = data.dNdX(i, 0); ///< exx
                B(1, i * 3 + 1) = data.dNdX(i, 1); ///< eyy
                B(2, i * 3 + 2) = data.dNdX(i, 2); ///< ezz

                B(3, i * 3 + 0) = data.dNdX(i, 1); ///< exy
                B(3, i * 3 + 1) = data.dNdX(i, 0); ///< exy

                B(4, i * 3 + 1) = data.dNdX(i, 2); ///< eyz
                B(4, i * 3 + 2) = data.dNdX(i, 1); ///< eyz

                B(5, i * 3 + 0) = data.dNdX(i, 2); ///< ezx
                B(5, i * 3 + 2) = data.dNdX(i, 0); ///< ezx
            }
        }
    }
    //=============================================================================
    void LinearElasticity::operator()(la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = U_.space();
        const auto &rho_ = constitutive_.rho();

        /// @todo The element matrices could be defined
        /// once here, and resized if nencessary, to
        /// reduce memory allocations

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &cell,
                        int cell_idx)
        {
            const auto element = V->element(cell.type);
            const auto int_rule = element->integration_rule();
            const int n_nodes = element->n_nodes();
            const int dim = element->dim();

            const auto elem_dof = V->cell_dof(cell_idx);
            const auto elem_pts = V->cell_dof_points(cell_idx);

            const int n_dof = n_nodes * dim;
            const int n_strain = constitutive_.n_strain();

            // Element displacement vector
            // la::DenseMatrix u(n_dof, 1);
            // Strain vector
            // la::DenseMatrix strain(n_strain, 1);
            // Stress vector
            // la::DenseMatrix stress(n_strain, 1);

            // Element strain-displacement matrix
            la::DenseMatrix B(n_strain, n_dof);

            // Element stress-strain matrix
            la::DenseMatrix D(n_strain, n_strain);

            // Element stiffness matrix
            la::DenseMatrix K(n_dof, n_dof);

            // Element force vector
            la::DenseMatrix F(n_dof, 1);

            for (int qpt_idx = 0; qpt_idx < int_rule->n_points(); qpt_idx++)
            {
                const real_t qwt = int_rule->weight(qpt_idx);
                const auto qpt = int_rule->point(qpt_idx);
                const auto data = element->transform(dim, qpt, elem_pts);
                const real_t Jwt = data.detJ * qwt;

                strain_displacement(data, B);

                constitutive_.tangent(cell_idx, cell.type, qpt, D);

                K += B.T() * D * B * Jwt;

                const real_t rho = rho_.value(cell_idx, cell.type, qpt);
                for (int i = 0; i < n_nodes; i++)
                {
                    // Inertial force
                    for (int dir = 0; dir < dim; dir++)
                    {
                        F(i * dim + dir, 0) += rho * g_[dir] * data.N(i, 0) * Jwt;
                    }

                    /// @todo Add user-defined force terms here
                }
            }

            lhs(elem_dof, elem_dof, K.values());
            if (rhs)
            {
                rhs(elem_dof, F.values());
            }
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }
}