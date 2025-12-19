#include "linear_elasticity.hpp"
#include <sfem/mesh/utils/loop_utils.hpp>
#include <sfem/mesh/utils/geo_utils.hpp>

namespace sfem::fem::solid_mechanics
{
    //=============================================================================
    LinearElasticity::LinearElasticity(FEField U, Strain &strain,
                                       LinearElasticIsotropic &constitutive,
                                       const std::array<real_t, 3> &g)
        : U_(U),
          strain_(strain),
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
    void LinearElasticity::operator()(la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = U_.space();
        const auto &rho_ = constitutive_.rho();
        const int dim = constitutive_.dim();
        const int n_strain = constitutive_.n_strain();

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
            const int n_dof = n_nodes * dim;

            const auto elem_dof = V->cell_dof(cell_idx);
            const auto elem_pts = V->cell_dof_points(cell_idx);

            // Element displacement vector
            // la::DenseMatrix u(n_dof, 1);

            // Strain vector
            la::DenseMatrix e(n_strain, 1);

            // Stress vector
            // la::DenseMatrix s(n_strain, 1);

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
                const auto data = element->transform(cell_idx, dim, qpt, elem_pts);
                const real_t Jwt = data.detJ * qwt;

                strain_.B_mat(data, B);

                constitutive_.tangent(data, e, D);

                K += B.T() * D * B * Jwt;

                const real_t rho = rho_.cell_value(cell_idx, qpt);
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