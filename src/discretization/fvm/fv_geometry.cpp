#include "fv_geometry.hpp"
#include "../fem/cg_space.hpp"
#include "../../mesh/utils/geo_utils.hpp"
#include "../../geo/utils.hpp"

namespace sfem::fvm
{
    FVGeometry::FVGeometry(std::shared_ptr<const mesh::Mesh> mesh)
        : mesh_(mesh)
    {
        // Quick access
        const auto topology = mesh_->topology();
        int dim = topology->dim();
        int n_cells = topology->n_entities(dim);
        int n_facets = topology->n_entities(dim - 1);

        // CG space for integration
        fem::CGSpace cg_space(mesh_, 1, {""});

        // Cell volumes and midpoints
        cell_volumes_.resize(n_cells, 0.0);
        cell_midpoints_.resize(n_cells);
        for (int i = 0; i < n_cells; i++)
        {
            auto cell_type = topology->entity(i, dim).type;
            auto cell_points = mesh_->entity_points(i, dim);
            auto element = cg_space.element(cell_type);
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
            auto facet_type = topology->entity(i, dim - 1).type;
            auto facet_points = mesh_->entity_points(i, dim - 1);
            auto element = cg_space.element(facet_type);
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
            std::array<geo::Vec3, 2> facet_cell_distance;
            for (int j = 0; j < 2; j++)
            {
                facet_cell_distance[j] = geo::compute_distance(cell_midpoints_[facet_adjacent_cells_[i][j]],
                                                               facet_midpoints_[i]);
                facet_cell_distances_[i][j] = facet_cell_distance[j].mag();
            }

            if (facet_adjacent_cells_[i][0] == facet_adjacent_cells_[i][1])
            {
                intercell_distances_[i] = 2 * facet_cell_distance[0];
            }
            else
            {
                intercell_distances_[i] = geo::Vec3(cell_midpoints_[facet_adjacent_cells_[i][0]],
                                                    cell_midpoints_[facet_adjacent_cells_[i][1]]);
            }
        }
    }
}