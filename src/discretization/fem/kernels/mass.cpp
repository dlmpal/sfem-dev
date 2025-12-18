#include "mass.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fem::kernels
{
    //=============================================================================
    MassND::MassND(FEField phi, Field &C)
        : phi_(phi),
          C_(C)
    {
    }
    //=============================================================================
    void MassND::operator()(la::MatSet lhs, la::VecSet)
    {
        // Quick access
        const auto V = phi_.space();
        const int n_comp = phi_.n_comp();

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

            // Element mass matrix
            la::DenseMatrix M(n_nodes * n_comp, n_nodes * n_comp);

            for (int qpt_idx = 0; qpt_idx < int_rule->n_points(); qpt_idx++)
            {
                const real_t qwt = int_rule->weight(qpt_idx);
                const auto qpt = int_rule->point(qpt_idx);
                const auto data = element->transform(cell_idx, dim, qpt, elem_pts);
                const real_t Jwt = data.detJ * qwt;

                const real_t C = C_.cell_value(cell_idx, qpt);
                for (int i = 0; i < n_nodes; i++)
                {
                    for (int j = 0; j < n_nodes; j++)
                    {
                        for (int k = 0; k < n_comp; k++)
                        {
                            M(i * n_comp + k, j * n_comp + k) += C * data.N(i, 0) * data.N(j, 0) * Jwt;
                        }
                    }
                }
            }
            lhs(elem_dof, elem_dof, M.values());
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }
}