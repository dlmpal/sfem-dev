#include "pressure_load.hpp"
#include <sfem/mesh/utils/geo_utils.hpp>
#include <sfem/mesh/utils/loop_utils.hpp>

namespace sfem::fem::solid_mechanics
{
    //=============================================================================
    PressureLoad::PressureLoad(FEField U, Field &P, const mesh::Region &region)
        : U_(U),
          P_(P),
          region_(region)
    {
    }
    //=============================================================================
    Field &PressureLoad::P()
    {
        return P_;
    }
    //=============================================================================
    const Field &PressureLoad::P() const
    {
        return P_;
    }
    //=============================================================================
    void PressureLoad::operator()(la::MatSet, la::VecSet rhs)
    {
        // Quick access
        const auto V = U_.space();
        const int dim = V->mesh()->pdim();

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &,
                        const mesh::Cell &facet,
                        int facet_idx)
        {
            const auto element = V->element(facet.type);
            const auto int_rule = element->integration_rule();
            const int n_nodes = element->n_nodes();
            const int n_dof = n_nodes * dim;

            const auto elem_dof = V->facet_dof(facet_idx);
            const auto elem_pts = V->facet_dof_points(facet_idx);
            const auto elem_normal = mesh::facet_normal(facet.type, elem_pts).normalize();

            la::DenseMatrix F(n_dof, 1);
            for (int qpt_idx = 0; qpt_idx < int_rule->n_points(); qpt_idx++)
            {
                const real_t qwt = int_rule->weight(qpt_idx);
                const auto qpt = int_rule->point(qpt_idx);
                const auto data = element->transform(facet_idx, dim, qpt, elem_pts);
                const real_t Jwt = data.detJ * qwt;

                const real_t pressure = P_.facet_value(facet_idx, qpt);
                for (int i = 0; i < n_nodes; i++)
                {
                    for (int dir = 0; dir < dim; dir++)
                    {
                        F(i * dim + dir, 0) += -pressure * data.N(i, 0) * elem_normal(dir) * Jwt;
                    }
                }
            }
            rhs(elem_dof, F.values());
        };
        mesh::utils::for_all_facets_region(*V->mesh(), work, region_);
    }
}