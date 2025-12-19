#include "post_utils.hpp"
#include <sfem/mesh/utils/loop_utils.hpp>

namespace sfem::fem
{
    //=============================================================================
    void cell_qpoint_average(const FESpace &V, ElementOperator op, FEField &F)
    {
        if (F.space()->order() > 0)
        {
            SFEM_ERROR("F should be a cell constant field\n");
        }

        // Quick access
        const auto mesh = V.mesh();
        const int dim = mesh->pdim();
        auto &F_values = F.dof_values();

        la::DenseMatrix f(F.n_comp(), 1);
        la::DenseMatrix fi(F.n_comp(), 1);

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &cell,
                        int cell_idx)
        {
            const auto element = V.element(cell.type);
            const auto int_rule = element->integration_rule();
            const auto elem_pts = V.cell_dof_points(cell_idx);

            f.set_all(0.0);
            fi.set_all(0.0);
            real_t vol = 0.0;
            for (int qpt_idx = 0; qpt_idx < int_rule->n_points(); qpt_idx++)
            {
                const real_t qwt = int_rule->weight(qpt_idx);
                const auto qpt = int_rule->point(qpt_idx);
                const auto data = element->transform(cell_idx, dim, qpt, elem_pts);
                const real_t Jwt = data.detJ * qwt;

                op(data, fi);
                f += fi * Jwt;
                vol += Jwt;
            }

            for (int comp_idx = 0; comp_idx < F.n_comp(); comp_idx++)
            {
                F_values(cell_idx, comp_idx) = f(comp_idx, 0) / std::abs(vol);
            }
        };
        mesh::utils::for_all_cells(*mesh, work);
        F_values.update_ghosts();
    }
}