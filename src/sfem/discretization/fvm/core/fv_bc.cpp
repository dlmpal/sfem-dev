#include "fv_bc.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FVBC::FVBC(const FVSpace &V, int n_comp)
        : n_comp_(n_comp)
    {
        // Initialize all boundary regions to zero-Neumann
        const auto mesh = V.mesh();
        int idx = 0;
        for (const auto &region : mesh->regions())
        {
            if (region.dim() < mesh->pdim())
            {
                const BCType region_type = BCType::zero_neumann;
                std::vector<int> region_facets;
                for (const auto &[facet, facet_idx] : mesh->region_facets(region.name()))
                {
                    region_facets.push_back(facet_idx);
                    bc_idx_[facet_idx] = idx++;
                }
                region_data_.insert({region.name(), {region_type, region_facets}});
            }
        }
        bc_data_.resize(idx * n_comp_);
    }
    //=============================================================================
    BCType FVBC::region_type(const std::string &region_name) const
    {
        return region_data_.at(region_name).first;
    }
    //=============================================================================
    std::vector<int> FVBC::region_facets(const std::string &region_name) const
    {
        return region_data_.at(region_name).second;
    }
    //=============================================================================
    void FVBC::set_region_bc(const std::string &region_name, BCType type, real_t value, int comp_idx)
    {
        set_region_bc(region_name, type, BCData{.c = value}, comp_idx);
    }
    //=============================================================================
    void FVBC::set_region_bc(const std::string &region_name, BCType type, BCData value, int comp_idx)
    {
        region_data_.at(region_name).first = type;
        for (const auto &facet_idx : region_data_.at(region_name).second)
        {
            const int bc_idx = bc_idx_[facet_idx];
            bc_data_[bc_idx * n_comp_ + comp_idx] = value;
        }
    }
    //=============================================================================
    real_t &FVBC::coeff(int facet_idx, int comp_idx)
    {
        return bc_data_.at(bc_idx_.at(facet_idx) * n_comp_ + comp_idx).a;
    }
    //=============================================================================
    real_t FVBC::coeff(int facet_idx, int comp_idx) const
    {
        return bc_data_.at(bc_idx_.at(facet_idx) * n_comp_ + comp_idx).a;
    }
    //=============================================================================
    real_t &FVBC::grad_coeff(int facet_idx, int comp_idx)
    {
        return bc_data_.at(bc_idx_.at(facet_idx) * n_comp_ + comp_idx).b;
    }
    //=============================================================================
    real_t FVBC::grad_coeff(int facet_idx, int comp_idx) const
    {
        return bc_data_.at(bc_idx_.at(facet_idx) * n_comp_ + comp_idx).b;
    }
    //=============================================================================
    real_t &FVBC::value(int facet_idx, int comp_idx)
    {
        return bc_data_.at(bc_idx_.at(facet_idx) * n_comp_ + comp_idx).c;
    }
    //=============================================================================
    real_t FVBC::value(int facet_idx, int comp_idx) const
    {
        return bc_data_.at(bc_idx_.at(facet_idx) * n_comp_ + comp_idx).c;
    }
}