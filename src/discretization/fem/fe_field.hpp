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
        /// reference point in a cell
        /// @param cell_idx Cell index
        /// @param cell_type Cell type
        /// @param pt Cell reference point
        /// @param comp_idx Component index
        /// @return Field value
        virtual real_t cell_value(int cell_idx,
                                  const std::array<real_t, 3> &pt,
                                  int comp_idx = 0) const = 0;

        /// @brief Evaluate the value of the field for a given
        /// reference point on a facet
        /// @param facet_idx Facet index
        /// @param pt Facet reference point
        /// @param comp_idx Component index
        /// @return Field value
        virtual real_t facet_value(int facet_idx,
                                   const std::array<real_t, 3> &pt,
                                   int comp_idx = 0) const = 0;

        /// @brief Evaluate the gradient of the field for a given
        /// reference point in a cell
        /// @param cell_idx Cell index
        /// @param cell_type Cell type
        /// @param pt Cell reference point
        /// @param comp_idx Component index
        /// @return Evaluated field gradient
        virtual geo::Vec3 cell_grad(int cell_idx,
                                    const std::array<real_t, 3> &pt,
                                    int comp_idx = 0) const = 0;

        /// @brief Evaluate the gradient of the field for a given
        /// reference point on a facet
        /// @param facet_idx Facet index
        /// @param pt Facet reference point
        /// @param comp_idx Component index
        /// @return Evaluated field gradient
        virtual geo::Vec3 facet_grad(int facet_idx,
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

        real_t cell_value(int cell_idx,
                          const std::array<real_t, 3> &pt,
                          int comp_idx = 0) const override;

        real_t facet_value(int facet_idx,
                           const std::array<real_t, 3> &pt,
                           int comp_idx = 0) const override;

        geo::Vec3 cell_grad(int cell_idx,
                            const std::array<real_t, 3> &pt,
                            int comp_idx = 0) const override;

        geo::Vec3 facet_grad(int facet_idx,
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

        real_t cell_value(int cell_idx,
                          const std::array<real_t, 3> &pt,
                          int comp_idx = 0) const override;

        real_t facet_value(int facet_idx,
                           const std::array<real_t, 3> &pt,
                           int comp_idx = 0) const override;

        geo::Vec3 cell_grad(int cell_idx,
                            const std::array<real_t, 3> &pt,
                            int comp_idx = 0) const override;

        geo::Vec3 facet_grad(int facet_idx,
                             const std::array<real_t, 3> &pt,
                             int comp_idx = 0) const override;

    private:
        /// @brief Finite element space
        std::shared_ptr<const FESpace> V_;

        /// @brief Mesh topology, obtained from the FE space's
        /// mesh. Stored here to reduce pointer dereferences
        std::shared_ptr<const mesh::Topology> topo_;

        /// @brief DoF values
        std::shared_ptr<la::Vector> dof_values_;
    };
}