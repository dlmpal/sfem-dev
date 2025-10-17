#include "fv_space.hpp"
#include "../fem/cg_space.hpp"
#include "../../mesh/utils/geo_utils.hpp"
#include "../../geo/utils.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FVSpace::FVSpace(std::shared_ptr<const mesh::Mesh> mesh)
        : mesh_(mesh)
    {
        // Quick access
        const int dim = mesh_->pdim();
        const auto topology = mesh_->topology();
        const int n_cells = topology->n_entities(dim);
        const int n_facets = topology->n_entities(dim - 1);
        const auto cell_to_facet = topology->connectivity(dim, dim - 1);
        const auto cell_index_map = topology->entity_index_map(dim);

        // DoF-to-DoF connectivity
        // The cell-to-cell conn. in topology is computed from the cell-to-node conn.
        // For the finite volume method, two cells are adjacent not when they share
        // a node, but a facet. Thus, the DoF-to-DoF connectivity is defined as follows
        connectivity_ = std::make_shared<graph::Connectivity>(cell_to_facet->primary_to_primary(1, true));

        // The cell index map in topology is not renumbered, i.e. the global indices of
        // the owned cells are not contiguous. Thus, the DoF index map is defined as
        // the renumbered cell index map
        index_map_ = std::make_shared<IndexMap>(topology->entity_index_map(dim)->renumber());

        // CG space for integration
        const fem::CGSpace cg_space(mesh_, 1);

        // Cell volumes and midpoints
        cell_volumes_.resize(n_cells, 0.0);
        cell_midpoints_.resize(n_cells);
        for (int i = 0; i < n_cells; i++)
        {
            const auto cell_type = topology->entity(i, dim).type;
            const auto cell_points = mesh_->entity_points(i, dim);
            const auto element = cg_space.element(cell_type);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                cell_volumes_[i] += element->transform(dim,
                                                       element->integration_rule()->point(nqpt),
                                                       cell_points)
                                        .detJ;
            }
            cell_midpoints_[i] = mesh::cell_midpoint(cell_points);
        }

        // Facet areas, midpoints, normal vectors and adjacent cells
        facet_areas_.resize(n_facets, 0.0);
        facet_midpoints_.resize(n_facets);
        facet_normals_.reserve(n_facets);
        facet_adjacent_cells_.resize(n_facets);
        for (int i = 0; i < n_facets; i++)
        {
            const auto facet_type = topology->entity(i, dim - 1).type;
            const auto facet_points = mesh_->entity_points(i, dim - 1);
            const auto element = cg_space.element(facet_type);
            for (int nqpt = 0; nqpt < element->integration_rule()->n_points(); nqpt++)
            {
                facet_areas_[i] += element->transform(dim,
                                                      element->integration_rule()->point(nqpt),
                                                      facet_points)
                                       .detJ;
            }
            facet_midpoints_[i] = mesh::cell_midpoint(facet_points);
            facet_normals_[i] = mesh::facet_normal(facet_type, facet_points).normalize();
            facet_adjacent_cells_[i] = topology->facet_adjacent_cells(i);
        }

        // Distances
        facet_cell_distances_.resize(n_facets);
        intercell_distances_.resize(n_facets);
        for (int i = 0; i < n_facets; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                facet_cell_distances_[i][j] = geo::compute_distance(cell_midpoints_[facet_adjacent_cells_[i][j]],
                                                                    facet_midpoints_[i]);
            }

            if (facet_adjacent_cells_[i][0] == facet_adjacent_cells_[i][1])
            {
                intercell_distances_[i] = 2 * geo::compute_distance(cell_midpoints_[facet_adjacent_cells_[i][0]],
                                                                    facet_midpoints_[i]);
            }
            else
            {
                intercell_distances_[i] = geo::Vec3(cell_midpoints_[facet_adjacent_cells_[i][0]],
                                                    cell_midpoints_[facet_adjacent_cells_[i][1]]);
            }
        }
    }
    //=============================================================================
    std::shared_ptr<const mesh::Mesh> FVSpace::mesh() const
    {
        return mesh_;
    }
    //=============================================================================
    std::shared_ptr<const graph::Connectivity> FVSpace::connectivity() const
    {
        return connectivity_;
    }
    //=============================================================================
    std::shared_ptr<const IndexMap> FVSpace::index_map() const
    {
        return index_map_;
    }
    //=============================================================================
    real_t FVSpace::cell_volume(int cell_idx) const
    {
        return cell_volumes_[cell_idx];
    }
    //=============================================================================
    real_t FVSpace::facet_area(int facet_idx) const
    {
        return facet_areas_[facet_idx];
    }
    //=============================================================================
    std::array<real_t, 3> FVSpace::cell_midpoint(int cell_idx) const
    {
        return cell_midpoints_[cell_idx];
    }
    //=============================================================================
    std::array<real_t, 3> FVSpace::facet_midpoint(int facet_idx) const
    {
        return facet_midpoints_[facet_idx];
    }
    //=============================================================================
    geo::Vec3 FVSpace::facet_normal(int facet_idx) const
    {
        return facet_normals_[facet_idx];
    }
    //=============================================================================
    geo::Vec3 FVSpace::intercell_distance(int facet_idx) const
    {
        return intercell_distances_[facet_idx];
    }
    //=============================================================================
    std::array<real_t, 2> FVSpace::facet_cell_distances(int facet_idx) const
    {
        return facet_cell_distances_[facet_idx];
    }
    //=============================================================================
    std::array<int, 2> FVSpace::facet_adjacent_cells(int facet_idx) const
    {
        return facet_adjacent_cells_[facet_idx];
    }
    //=============================================================================
    bool FVSpace::is_boundary(int facet_idx) const
    {
        const auto &cells = facet_adjacent_cells_[facet_idx];
        if (cells[0] == cells[1])
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    //=============================================================================
    real_t FVSpace::compute_facet_value(int facet_idx, real_t value1, real_t value2, bool harmonic) const
    {
        const real_t d = intercell_distances_[facet_idx].mag();
        const auto &[d1, d2] = facet_cell_distances_[facet_idx];
        const real_t w = d2 / d;
        if (harmonic)
        {
            return 1 / (w / value1 + (1 - w) / value2);
        }
        else
        {
            return w * value1 + (1 - w) * value2;
        }
    }
}