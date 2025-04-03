#include "fe_space.hpp"
#include <set>

namespace sfem::fem
{
    //=============================================================================
    FESpace::FESpace(std::shared_ptr<const mesh::Mesh> mesh, int order,
                     const std::vector<std::string> &components,
                     const std::string &name)
        : mesh_(mesh),
          order_(order),
          components_(components),
          name_(name)
    {
        if (order_ <= 0)
        {
            SFEM_ERROR(std::format("Cannot create FESpace with order {} (<=0)\n", order_));
        }
        if (components_.size() <= 0)
        {
            SFEM_ERROR(std::format("Cannot create FESpace with {} (<=0) components\n",
                                   components_.size()));
        }
    }
    //=============================================================================
    FESpace::~FESpace()
    {
    }
    //=============================================================================
    std::shared_ptr<const mesh::Mesh> FESpace::mesh() const
    {
        return mesh_;
    }
    //=============================================================================
    int FESpace::order() const
    {
        return order_;
    }
    //=============================================================================
    std::vector<std::string> FESpace::components() const
    {
        return components_;
    }
    //=============================================================================
    std::string FESpace::name() const
    {
        return name_;
    }
    //=============================================================================
    const std::array<std::shared_ptr<graph::Connectivity>, 2> &FESpace::connectivity() const
    {
        return connectivity_;
    }
    //=============================================================================
    std::shared_ptr<const IndexMap> FESpace::index_map() const
    {
        return index_map_;
    }
    //=============================================================================
    const FECollection &FESpace::fe_collection() const
    {
        return fe_collection_;
    }
    //=============================================================================
    int FESpace::n_comp() const
    {
        return static_cast<int>(components_.size());
    }
    //=============================================================================
    int FESpace::comp_idx(const std::string &component) const
    {
        auto it = std::find(components_.cbegin(),
                            components_.end(),
                            component);
        if (it == components_.end())
        {
            auto msg = std::format("Component {} not found in: [ ", component);
            for (const auto &comp : components_)
            {
                msg += std::format("{} ", comp);
            }
            msg += "]\n";
            log_msg(msg, LogLevel::warning);
            return -1;
        }
        else
        {
            return static_cast<int>(std::distance(components_.cbegin(), it));
        }
    }
    //=============================================================================
    std::shared_ptr<FiniteElement> FESpace::element(mesh::CellType cell_type) const
    {
        auto element = fe_collection_[static_cast<int>(cell_type)];
        if (element == nullptr)
        {
            SFEM_ERROR(std::format("Cell type {} is not supported for {} space\n",
                                   mesh::cell_type_str(cell_type),
                                   name_));
        }
        return element;
    }
    //=============================================================================
    std::span<const int> FESpace::cell_dof(int cell_idx) const
    {
        return connectivity_[0]->links(cell_idx);
    }
    //=============================================================================
    std::vector<int> FESpace::facet_dof(int facet_idx) const
    {
        // Quick access
        const auto topology = mesh_->topology();
        int dim = topology->dim();

        // Owner cell index and facet relative index (in the owner cell)
        auto owner_cell_idx = topology->entity_owner(facet_idx, dim - 1);
        auto facet_rel_idx = topology->entity_rel_idx(owner_cell_idx, dim,
                                                      facet_idx, dim - 1);

        // Owner cell and facet types
        auto owner_cell_type = topology->entity(owner_cell_idx, dim).type;
        auto facet_type = topology->entity(facet_idx, dim - 1).type;

        // Owner cell and facet DoF
        auto owner_cell_dof = connectivity_[0]->links(owner_cell_idx);
        std::vector<int> facet_dof(dof::cell_num_dof(facet_type, order_));
        if (dim == 3)
        {
            auto face_ordering = mesh::cell_face_ordering(owner_cell_type, facet_rel_idx);
            for (int i = 0; i < mesh::cell_num_nodes(facet_type); i++)
            {
                facet_dof[i] = owner_cell_dof[face_ordering[i]];
            }
        }
        else if (dim == 2)
        {
            auto edge_ordering = mesh::cell_edge_ordering(owner_cell_type, facet_rel_idx);
            for (int i = 0; i < 2; i++)
            {
                facet_dof[i] = owner_cell_dof[edge_ordering[i]];
            }
            int offset = facet_rel_idx * dof::cell_num_internal_dof(facet_type, order_) +
                         mesh::cell_num_nodes(owner_cell_type);
            for (int i = 0; i < dof::cell_num_internal_dof(facet_type, order_); i++)
            {
                facet_dof[i + 2] = owner_cell_dof[offset + i];
            }
        }
        else
        {
            facet_dof[0] = owner_cell_dof[facet_rel_idx];
        }
        return facet_dof;
    }
    //=============================================================================
    std::vector<int> FESpace::boundary_dof(const std::string &region_name) const
    {
        // Get all DoFs belonging to the boundary region
        std::set<int> boundary_dof_;
        for (const auto &[_, facet_idx] : mesh_->region_facets(region_name))
        {
            for (int dof : facet_dof(facet_idx))
            {
                boundary_dof_.insert(dof);
            }
        }

        // Copy the DoFs into a vector
        std::vector<int> boundary_dof(boundary_dof_.size());
        std::copy(boundary_dof_.cbegin(),
                  boundary_dof_.cend(),
                  boundary_dof.begin());

        return boundary_dof;
    }
    //=============================================================================
    std::vector<std::array<real_t, 3>> FESpace::cell_dof_points(int cell_idx) const
    {
        auto cell = mesh_->topology()->entity(cell_idx, mesh_->topology()->dim());
        auto points = mesh_->entity_points(cell_idx, mesh_->topology()->dim());
        dof::compute_cell_dof_points(cell.type, order_, points);
        return points;
    }
    //=============================================================================
    std::vector<std::array<real_t, 3>> FESpace::facet_dof_points(int facet_idx) const
    {
        auto facet = mesh_->topology()->entity(facet_idx, mesh_->topology()->dim() - 1);
        auto points = mesh_->entity_points(facet_idx, mesh_->topology()->dim() - 1);
        dof::compute_cell_dof_points(facet.type, order_, points);
        return points;
    }
    //=============================================================================
    std::vector<std::array<real_t, 3>> FESpace::dof_points() const
    {
        std::vector<std::array<real_t, 3>> points(connectivity_[0]->n_secondary());
        for (int i = 0; i < connectivity_[0]->n_primary(); i++)
        {
            auto dof = cell_dof(i);
            auto dof_points = cell_dof_points(i);
            for (int j = 0; j < connectivity_[0]->n_links(i); j++)
            {
                points[dof[j]] = dof_points[j];
            }
        }
        return points;
    }
}