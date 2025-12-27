#include "mesh.hpp"
#include <sfem/base/error.hpp>
#include <algorithm>
#include <ranges>
#include <format>

namespace sfem::mesh
{
    //=============================================================================
    Mesh::Mesh(std::shared_ptr<Topology> topology,
               std::vector<std::array<real_t, 3>> &&points,
               std::vector<Region> &&regions)
        : topology_(topology),
          points_(std::move(points)),
          regions_(std::move(regions))
    {
        SFEM_CHECK_SIZES(topology_->n_entities(0), points_.size());

        // Compute the mesh's physical dimension as the max dim. of the regions
        // Require that the topological dimension is equal, or less than the
        // physical
        dim_ = std::ranges::max(regions_,
                                [](Region &lhs, Region &rhs)
                                { return lhs.dim() < rhs.dim(); })
                   .dim();
        if (topology_->dim() > dim_)
        {
            SFEM_ERROR(std::format("Topological dimension ({}) is greater than physical dimension ({})\n",
                                   topology_->dim(),
                                   dim_));
        }
    }
    //=============================================================================
    std::shared_ptr<Topology> Mesh::topology()
    {
        return topology_;
    }
    //=============================================================================
    std::shared_ptr<const Topology> Mesh::topology() const
    {
        return topology_;
    }
    //=============================================================================
    std::vector<std::array<real_t, 3>> &Mesh::points()
    {
        return points_;
    }
    //=============================================================================
    const std::vector<std::array<real_t, 3>> &Mesh::points() const
    {
        return points_;
    }
    //=============================================================================
    const std::vector<Region> &Mesh::regions() const
    {
        return regions_;
    }
    //=============================================================================
    int Mesh::pdim() const
    {
        return dim_;
    }
    //=============================================================================
    int Mesh::tdim() const
    {
        return topology_->dim();
    }
    //=============================================================================
    std::vector<std::array<real_t, 3>>
    Mesh::entity_points(int entity_idx, int dim) const
    {
        auto entity_nodes = topology_->adjacent_entities(entity_idx, dim, 0);
        std::vector<std::array<real_t, 3>> entity_points(entity_nodes.size());
        for (std::size_t i = 0; i < entity_points.size(); i++)
        {
            entity_points[i] = points_[entity_nodes[i]];
        }
        return entity_points;
    }
    //=============================================================================
    Region Mesh::get_region_by_name(const std::string &region_name) const
    {
        const auto it = std::find_if(regions_.cbegin(),
                                     regions_.cend(),
                                     [&region_name](const Region &region)
                                     { return region.name() == region_name; });
        if (it == regions_.end())
        {
            SFEM_ERROR(std::format("Invalid region name: {} \n", region_name));
        }
        return *it;
    }
    //=============================================================================
    Region Mesh::get_region_by_tag(int region_tag) const
    {
        const auto it = std::find_if(regions_.cbegin(),
                                     regions_.cend(),
                                     [region_tag](const Region &region)
                                     {
                                         return region.tag() == region_tag;
                                     });
        if (it == regions_.end())
        {
            SFEM_ERROR(std::format("Invalid region tag: {} \n", region_tag));
        }
        return *it;
    }
    //=============================================================================
    std::vector<std::pair<Cell, int>>
    Mesh::region_cells(const std::string &region_name) const
    {
        auto region = get_region_by_name(region_name);
        std::vector<std::pair<Cell, int>> region_cells;
        for (int i = 0; i < topology_->n_entities(topology_->dim()); i++)
        {
            if (auto cell = topology_->entity(i, topology_->dim()); cell.tag == region.tag())
            {
                region_cells.emplace_back(std::make_pair(cell, i));
            }
        }
        return region_cells;
    }
    //=============================================================================
    std::vector<std::pair<Cell, int>>
    Mesh::region_facets(const std::string &region_name) const
    {
        auto region = get_region_by_name(region_name);
        std::vector<std::pair<Cell, int>> region_facets;
        for (int i = 0; i < topology_->n_entities(topology_->dim() - 1); i++)
        {
            if (auto facet = topology_->entity(i, topology_->dim() - 1); facet.tag == region.tag())
            {
                region_facets.emplace_back(std::make_pair(facet, i));
            }
        }
        return region_facets;
    }
}