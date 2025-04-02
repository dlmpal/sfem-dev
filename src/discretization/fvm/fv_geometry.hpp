#pragma once

#include "../../mesh/mesh.hpp"

namespace sfem::fvm
{
    class FVGeometry
    {
    public:
        FVGeometry(std::shared_ptr<const mesh::Mesh> mesh);

        std::shared_ptr<const mesh::Mesh> mesh_;

        std::vector<real_t> cell_volumes_;
        std::vector<real_t> facet_areas_;

        std::vector<std::array<real_t, 3>> cell_midpoints_;
        std::vector<std::array<real_t, 3>> facet_midpoints_;

        std::vector<geo::Vec3> facet_normals_;
        std::vector<geo::Vec3> intercell_distances_;
        std::vector<std::array<real_t, 2>> facet_cell_distances_;
        std::vector<std::array<int, 2>> facet_adjacent_cells_;
    };
}