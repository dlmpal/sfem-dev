#pragma once

#include "fe_space.hpp"
#include "../../la/native/vector.hpp"

namespace sfem::fem
{
    class Field
    {
    public:
        Field(const std::vector<std::string> &components);
        virtual ~Field() = default;

        int n_comp() const;
        std::vector<std::string> components() const;

        /// @brief Evaluate the value of the field for a given
        /// reference point inside a cell
        /// @param cell_idx Cell index
        /// @param cell_type Cell type
        /// @param pt Reference point
        /// @param comp_idx Component index
        /// @return Evaluated field value
        virtual real_t value(int cell_idx,
                             mesh::CellType cell_type,
                             const std::array<real_t, 3> &pt,
                             int comp_idx = 0) const = 0;

        /// @brief Evaluate the gradient of the field for a given
        /// reference point inside a cell
        /// @param cell_idx Cell index
        /// @param cell_type Cell type
        /// @param pt Reference point
        /// @param comp_idx Component index
        /// @return Evaluated field gradient
        virtual geo::Vec3 grad(int cell_idx,
                               mesh::CellType cell_type,
                               const std::array<real_t, 3> &pt,
                               int comp_idx = 0) const = 0;

    protected:
        std::vector<std::string> components_;
    };

    class ConstantField : public Field
    {
    public:
        ConstantField(const std::vector<std::string> &components,
                      const std::vector<real_t> &value);

        ConstantField(const std::string &name, real_t value);

        real_t value(int cell_idx,
                     mesh::CellType cell_type,
                     const std::array<real_t, 3> &pt,
                     int comp_idx = 0) const override;

        geo::Vec3 grad(int cell_idx,
                       mesh::CellType cell_type,
                       const std::array<real_t, 3> &pt,
                       int comp_idx = 0) const override;

    private:
        std::vector<real_t> value_;
    };

    class FEField : public Field
    {
    public:
        FEField(std::shared_ptr<const FESpace> V,
                const std::vector<std::string> &components);

        std::shared_ptr<const FESpace> space() const;

        la::Vector &dof_values();
        const la::Vector &dof_values() const;

        void cell_values(int cell_idx,
                         std::span<real_t> values) const;

        void cell_values(int cell_idx, int comp_idx,
                         std::span<real_t> values) const;

        real_t value(int cell_idx,
                     mesh::CellType cell_type,
                     const std::array<real_t, 3> &pt,
                     int comp_idx = 0) const override;

        geo::Vec3 grad(int cell_idx,
                       mesh::CellType cell_type,
                       const std::array<real_t, 3> &pt,
                       int comp_idx = 0) const override;

    private:
        std::shared_ptr<const FESpace> V_;

        std::shared_ptr<la::Vector> dof_values_;
    };
}