#pragma once

#include <sfem/discretization/fem/core/fe_space.hpp>

namespace sfem::fem
{
    /// @brief This class enables efficient and convenient handling
    /// of Dirichlet boundary conditions for a field defined on a
    // finite-element space
    class DirichletBC
    {
    public:
        /// @brief Create a DirichletBC
        /// @param phi Finite-element field
        DirichletBC(std::shared_ptr<const FESpace> V, int n_comp = 1);

        // Avoid copies
        DirichletBC(const DirichletBC &) = delete;
        DirichletBC &operator=(const DirichletBC &) = delete;

        // Move constructor and assignment
        DirichletBC(DirichletBC &&) = default;
        DirichletBC &operator=(DirichletBC &&) = default;

        /// @brief Set the value of a specific component on a boundary region
        void set_value(const std::string &region_name, real_t value, int comp_idx = 0);

        /// @brief Set the values of a specific component on a boundary region
        void set_values(const std::string &region_name,
                        std::span<const real_t> values,
                        int comp_idx = 0);

        /// @brief Get all specified DoFs and their values
        std::pair<std::vector<int>, std::vector<real_t>>
        get_dofs_values() const;

        /// @brief Reset/clears all specified values
        void reset_values();

    private:
        std::shared_ptr<const FESpace> V_;

        int n_comp_;

        /// @brief The DoF belonging to boundary regions.
        /// They are obtained from the field's FESpace only when a
        /// B.C. is specified on a boundary region for the first time.
        /// Thus, re-specifying the values on a boundary is relatively cheap
        std::unordered_map<std::string, std::vector<int>> boundary_dof_;

        std::unordered_map<int, real_t> data_;
    };
}
