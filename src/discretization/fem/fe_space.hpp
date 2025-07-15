#pragma once

#include "elements/fe.hpp"
#include "../../mesh/mesh.hpp"

namespace sfem::fem
{
    /// @brief A collection of finite elements assigns an element to each suitable cell type
    using FECollection = std::array<std::shared_ptr<FiniteElement>,
                                    static_cast<std::size_t>(mesh::CellType::n_cell_types)>;

    /// @brief Finite element space ABC.
    /// Each subclass should define the cell-to-DoF connectivity,
    /// the DoF index map and the finite element collection in its
    /// constructor
    class FESpace
    {
    public:
        /// @brief Create a finite element space
        /// @param mesh Mesh
        /// @param order Polynomial order (or degree)
        /// @param name The name of the space
        FESpace(std::shared_ptr<const mesh::Mesh> mesh,
                int order,
                const std::string &name);

        // Destructor
        virtual ~FESpace() = 0;

        /// @brief Get the underlying mesh
        std::shared_ptr<const mesh::Mesh> mesh() const;

        /// @brief Get the order of the finite element space
        int order() const;

        /// @brief Get the name of the finite element space
        std::string name() const;

        /// @brief Get the cell-to-DoF and DoF-to-DoF connectivities
        const std::array<std::shared_ptr<graph::Connectivity>, 2> &connectivity() const;

        /// @brief Get the DoF index map
        std::shared_ptr<const IndexMap> index_map() const;

        /// @brief Get the finite element collection
        const FECollection &fe_collection() const;

        /// @brief Get the finite element for the corresponding cell type
        /// @note Raises an error if no element exists for the given cell type
        std::shared_ptr<FiniteElement> element(mesh::CellType cell_type) const;

        /// @brief Get the DoF for the cell
        std::span<const int> cell_dof(int cell_idx) const;

        /// @brief Get the DoF for the facet
        std::vector<int> facet_dof(int facet_idx) const;

        /// @brief Get the DoF belonging to a boundary region
        std::vector<int> boundary_dof(const std::string &region_name) const;

        /// @brief Get the coordinates of the DoF for a cell
        std::vector<std::array<real_t, 3>> cell_dof_points(int cell_idx) const;

        /// @brief Get the coordinates of the DoF for a facet
        std::vector<std::array<real_t, 3>> facet_dof_points(int facet_idx) const;

        /// @brief Get the coordinates of all DoF
        std::vector<std::array<real_t, 3>> dof_points() const;

    protected:
        /// @brief The mesh
        std::shared_ptr<const mesh::Mesh> mesh_;

        /// @brief The order of the space (i.e. the polynomial degree)
        int order_;

        /// @brief The name of the space
        std::string name_;

        /// @brief The cell-to-DoF and DoF-to-DoF connectivity.
        /// The former is required by this class to perform various operations,
        /// while the latter is used to contruct sparse matrices.
        /// @note The DoF-to-DoF connectivity should include diagonal terms
        /// (i.e. each DoF should be connected to itself)
        std::array<std::shared_ptr<graph::Connectivity>, 2> connectivity_;

        /// @brief The DoF index map
        std::shared_ptr<IndexMap> index_map_;

        /// @brief The finite element collection,
        /// i.e. a finite element corresponding to each suitable cell type
        FECollection fe_collection_;
    };
}