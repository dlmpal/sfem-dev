#pragma once

#include "fv_space.hpp"

namespace sfem::fvm
{
    enum class BCType
    {
        zero_neumann,
        neumann,
        dirichlet,
        robin
    };

    /// @brief This struct represents different B.C. types
    /// as follows:
    /// 1) Dirichlet: u = c
    /// 2) Neumann: du/dn = c
    /// 3) Robin: a * du/dn + b * u = c
    struct BCData
    {
        real_t a = 0.0;
        real_t b = 0.0;
        real_t c = 0.0;
    };

    class FVBC
    {
    public:
        FVBC(const FVSpace &V, int n_comp);

        BCType region_type(const std::string &region_name) const;

        std::vector<int> region_facets(const std::string &region_name) const;

        void set_region_bc(const std::string &region_name, BCType type, real_t value, int comp_idx = 0);
        void set_region_bc(const std::string &region_name, BCType type, BCData value, int comp_idx = 0);

        real_t &facet_value(int facet_idx, int comp_idx = 0);
        real_t facet_value(int facet_idx, int comp_idx = 0) const;

        BCData &facet_data(int facet_idx, int comp_idx = 0);
        BCData facet_data(int facet_idx, int comp_idx = 0) const;

    private:
        int n_comp_;

        /// @todo Move this to FVSpace, store a pointer to the FVSpace here
        std::unordered_map<std::string, std::pair<BCType, std::vector<int>>> region_data_;

        std::unordered_map<int, int> bc_idx_;

        std::vector<BCData> bc_data_;
    };
}