#include "diffusion.hpp"
#include "../../../mesh/utils/loop_utils.hpp"

namespace sfem::fem::kernels
{
    //=============================================================================
    Diffusion::Diffusion(FEField phi, Field &D)
        : phi_(phi),
          D_(D)
    {
    }
    //=============================================================================
    FEField Diffusion::field() const
    {
        return phi_;
    }
    //=============================================================================
    Field &Diffusion::D()
    {
        return D_;
    }
    //=============================================================================
    const Field &Diffusion::D() const
    {
        return D_;
    }
    //=============================================================================
    void Diffusion::operator()(la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = phi_.space();

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

            la::DenseMatrix K(element->n_nodes(), element->n_nodes());

            for (int qpt_idx = 0; qpt_idx < int_rule->n_points(); qpt_idx++)
            {
                const real_t qwt = int_rule->weight(qpt_idx);
                const auto qpt = int_rule->point(qpt_idx);
                const auto data = element->transform(dim, qpt, elem_pts);
                const real_t Jwt = data.detJ * qwt;

                const real_t D = D_.value(cell_idx, cell.type, qpt);
                for (int i = 0; i < n_nodes; i++)
                {
                    for (int j = 0; j < n_nodes; j++)
                    {
                        real_t aij = 0.0;
                        for (int dir = 0; dir < dim; dir++)
                        {
                            aij += data.dNdX(i, dir) * data.dNdX(j, dir);
                        }
                        K(i, j) += D * aij * data.detJ * Jwt;
                    }
                }
            }
            lhs(elem_dof, elem_dof, K.values());
        };
        mesh::utils::for_all_cells(*V->mesh(), work);
    }
}