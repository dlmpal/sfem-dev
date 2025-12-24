#include "laplacian.hpp"
#include <sfem/discretization/fvm/core/fv_bc.hpp>
#include <sfem/la/native/setval_utils.hpp>
#include <sfem/mesh/utils/loop_utils.hpp>

namespace sfem::fvm
{
    //=============================================================================
    Laplacian::Laplacian(FVField phi, IField &D)
        : phi_(phi),
          D_(D)
    {
        if (phi_.n_comp() > 1)
        {
            SFEM_ERROR(std::format("Only scalar functions are supported, (n_comp={} > 1)\n", phi_.n_comp()));
        }
    }
    //=============================================================================
    FVField &Laplacian::field()
    {
        return phi_;
    }
    //=============================================================================
    const FVField &Laplacian::field() const
    {
        return phi_;
    }
    //=============================================================================
    IField &Laplacian::D()
    {
        return D_;
    }
    //=============================================================================
    const IField &Laplacian::D() const
    {
        return D_;
    }
    //=============================================================================
    void Laplacian::operator()(la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = phi_.space();
        const auto bc = phi_.boundary_condition();

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &region,
                        const mesh::Cell &,
                        int facet_idx)
        {
            // Cells adjacent to the facet
            const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
            const auto &[owner, neighbour] = adjacent_cells;

            // Facet area vector and intercell distance
            const geo::Vec3 Sf = V->facet_area_vec(facet_idx);
            const geo::Vec3 dPN = V->facet_intercell_distance(facet_idx);

            // Boundary facets
            if (owner == neighbour)
            {
                const std::array<int, 1> idx = {owner};
                std::array<real_t, 1> lhs_value{};
                std::array<real_t, 1> rhs_value{};

                // Facet diffusivity coefficient
                const real_t Df = D_.cell_value(owner);

                if (bc.region_type(region.name()) == BCType::dirichlet)
                {
                    lhs_value[0] = 2.0 * Df * Sf.mag() / dPN.mag();
                    rhs_value[0] = lhs_value[0] * bc.facet_value(facet_idx);
                }
                else if (bc.region_type(region.name()) == BCType::neumann)
                {
                    lhs_value[0] = 0.0;
                    rhs_value[0] = Df * Sf.mag() * bc.facet_value(facet_idx);
                }
                else if (bc.region_type(region.name()) == BCType::robin)
                {
                    const auto [a, b, c] = bc.facet_data(facet_idx);
                    const real_t h_inf = b / a;
                    const real_t phi_inf = c / a;
                    const real_t d = dPN.mag();
                    lhs_value[0] = (h_inf * Df / d) /
                                   (h_inf + Df / d) * Sf.mag();
                    rhs_value[0] = lhs_value[0] * phi_inf;
                }

                lhs(idx, idx, lhs_value);
                rhs(idx, rhs_value);
            }
            // Internal facets
            else
            {
                // Decompose the area vector
                const auto [delta, kappa] = V->decompose_area_vec(facet_idx);

                // Facet diffusivity coefficient
                const real_t Df = D_.facet_value(facet_idx);

                // LHS - orthogonal contribution
                {
                    const real_t value = Df * delta.mag() / dPN.mag();
                    std::array<real_t, 4> lhs_values = {value, -value,
                                                        -value, value};
                    lhs(adjacent_cells, adjacent_cells, lhs_values);
                }

                // RHS - non-orthogonal correction
                {
                    const real_t value = Df * geo::inner(phi_.facet_grad(facet_idx), kappa);
                    std::array<real_t, 2> rhs_values{value, -value};
                    rhs(adjacent_cells, rhs_values);
                }
            }
        };
        mesh::utils::for_all_facets(*V->mesh(), work);
    }
}