#pragma once

#include "../mesh/mesh.hpp"
#include "../la/native/vector.hpp"

namespace sfem
{
    class Function
    {
    public:
        Function(int n_comp);

        // Disable copying
        Function(const Function &) = delete;
        Function &operator=(const Function &) = delete;

        // Move constructor and assignment
        Function(Function &&) = default;
        Function &operator=(Function &&) = default;

        /// @brief Get the number of components
        int n_comp() const;

        virtual real_t &operator()(int idx, int comp) = 0;
        virtual real_t operator()(int idx, int comp) const = 0;

    protected:
        /// @brief Number of components
        int n_comp_;
    };

    class ConstantFunction : public Function
    {
    public:
        ConstantFunction(real_t value);

        real_t &value();
        real_t value() const;

        real_t &operator()(int idx, int comp);

        real_t operator()(int idx, int comp) const;

    private:
        /// @brief Constant value
        real_t value_;
    };

    class MeshFunction : public Function
    {
    public:
        MeshFunction(std::shared_ptr<mesh::Mesh> mesh, int dim, int n_comp);

        std::shared_ptr<const mesh::Mesh> mesh() const;

        la::Vector &values();
        const la::Vector &values() const;

        real_t &operator()(int idx, int comp);
        real_t operator()(int idx, int comp) const;

    private:
        /// @brief Mesh
        std::shared_ptr<mesh::Mesh> mesh_;

        /// @brief Values
        la::Vector values_;
    };
}