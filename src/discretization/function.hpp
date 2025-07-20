#pragma once

#include "../mesh/mesh.hpp"
#include "../la/native/vector.hpp"

namespace sfem
{
    class Function : public la::Vector
    {
    public:
        /// @param components The name of each component
        Function(std::shared_ptr<const IndexMap> index_map,
                 const std::vector<std::string> &components);

        /// @brief Get the name of each component
        std::vector<std::string> components() const;

        /// @brief Get the number of components
        int n_comp() const;

        /// @brief Get the index of a specific component
        /// @note Components are numbered according to their
        /// names were provided. For example, components {"u", "v"}
        /// have an index of 0 and 1 respectively
        /// @note If the component is not found, returns -1
        int comp_idx(const std::string &component) const;

    protected:
        /// @brief Names of function components
        std::vector<std::string> components_;
    };

    class MeshFunction : public Function
    {
    public:
        MeshFunction(std::shared_ptr<const mesh::Mesh> mesh, int dim,
                     const std::vector<std::string> &components);

        std::shared_ptr<const mesh::Mesh> mesh() const;

    private:
        /// @brief Mesh
        std::shared_ptr<const mesh::Mesh> mesh_;
    };
}