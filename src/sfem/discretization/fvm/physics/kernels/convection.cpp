#include "convection.hpp"
#include <sfem/discretization/fvm/core/fv_bc.hpp>
#include <sfem/la/native/setval_utils.hpp>
#include <sfem/mesh/utils/loop_utils.hpp>

namespace sfem::fvm
{
    //=============================================================================
    Convection::Convection(FVField phi, std::vector<real_t> &flux)
        : phi_(phi),
          flux_(flux)
    {
    }
    //=============================================================================
    FVField &Convection::field()
    {
        return phi_;
    }
    //=============================================================================
    const FVField &Convection::field() const
    {
        return phi_;
    }
    //=============================================================================
    std::vector<real_t> &Convection::flux()
    {
        return flux_;
    }
    //=============================================================================
    const std::vector<real_t> &Convection::flux() const
    {
        return flux_;
    }
    //=============================================================================
    void Convection::operator()(la::MatSet lhs, la::VecSet rhs)
    {
        // Quick access
        const auto V = phi_.space();
        const auto &bc = phi_.boundary_condition();

        auto work = [&](const mesh::Mesh &,
                        const mesh::Region &region,
                        const mesh::Cell &,
                        int facet_idx)
        {
            // Cells adjacent to the facet
            const auto adjacent_cells = V->facet_adjacent_cells(facet_idx);
            const auto &[owner, neighbour] = adjacent_cells;

            // Facet flux
            const real_t Ff = flux_[facet_idx];

            // Boundary facets
            if (owner == neighbour)
            {
                const std::array<int, 1> idx = {owner};
                std::array<real_t, 1> lhs_value{};
                std::array<real_t, 1> rhs_value{};

                if (bc.region_type(region.name()) == BCType::dirichlet)
                {
                    if (Ff >= 0)
                    {
                        lhs_value[0] = Ff;
                    }
                    else
                    {
                        const real_t phi_facet = bc.facet_value(facet_idx);
                        rhs_value[0] = -Ff * phi_facet;
                    }
                }
                else if (bc.region_type(region.name()) == BCType::neumann)
                {
                    if (Ff > 0)
                    {
                        lhs_value[0] = Ff;
                    }
                    else
                    {
                        const real_t grad_facet = bc.facet_value(facet_idx);
                        const real_t dfP = V->facet_cell_distances(facet_idx)[0];
                        lhs_value[0] = Ff;
                        rhs_value[0] = -Ff * dfP * grad_facet;
                    }
                }
                else if (bc.region_type(region.name()) == BCType::robin)
                {
                    /// @todo
                }
                else // Zero-Neumann
                {
                    lhs_value[0] = Ff;
                }

                lhs(idx, idx, lhs_value);
                rhs(idx, rhs_value);
            }
            // Internal facets
            else
            {
                // Upwind differencing
                const real_t w = Ff > 0 ? 1.0 : 0.0;
                const std::array<real_t, 4> lhs_values = {w * Ff, (1 - w) * Ff,
                                                          -w * Ff, -(1 - w) * Ff};
                lhs(adjacent_cells, adjacent_cells, lhs_values);
            }
        };
        mesh::utils::for_all_facets(*V->mesh(), work);
    }
}