#pragma once

#include "../../parallel/index_map.hpp"
#include "../../base/config.hpp"
#include <memory>

namespace sfem::la
{
    /// @brief MPI-parallel vector
    class Vector
    {
    public:
        /// @brief Create a vector
        /// @param index_map Index map
        /// @param block_size Block size (i.e. no. components per index)
        /// @param values Local values (owned + ghost)
        Vector(std::shared_ptr<const IndexMap> index_map, int block_size, std::vector<real_t> &&values);

        /// @brief Create a vector
        /// @param index_map Index map
        /// @param block_size Block size (i.e. no. components per index)
        /// @param value Uniform value
        Vector(std::shared_ptr<const IndexMap> im, int block_size, real_t value = 0.0);

        // Avoid uninentional copying
        Vector(const Vector &) = delete;
        Vector &operator=(const Vector &) = delete;

        // Move constructor and assignment operator
        Vector(Vector &&) = default;
        Vector &operator=(Vector &&) = default;

        /// @brief Get the vector's index map
        std::shared_ptr<const IndexMap> index_map() const;

        /// @brief Get the vector's block size
        int block_size() const;

        /// @brief Get the vector's values
        std::vector<real_t> &values();

        /// @brief Get the vector's values (const version)
        const std::vector<real_t> &values() const;

        /// @brief Get the vector's local size
        int n_local() const;

        /// @brief Get the vector's global size
        int n_global() const;

        /// @brief Get the value for a given (local) index and component
        real_t &operator()(int idx, int comp = 0);

        /// @brief Get the value for a given (local) index and component  (const version)
        real_t operator()(int idx, int comp = 0) const;

        /// @brief Set all vector values to a uniform value
        void set_all(real_t value);

        /// @brief Set vector values for a given set of (local) indices
        /// @param idxs Local indices
        /// @param values Values
        /// @param insert Whether to insert, or add the values
        void set_values(std::span<const int> idxs,
                        std::span<const real_t> values,
                        bool insert = false);

        /// @brief Assemble the vector, i.e. synchronize the values of ghost indices
        void assemble();

    protected:
        /// @brief Index map
        std::shared_ptr<const IndexMap> index_map_;

        /// @brief Block size
        int block_size_;

        /// @brief Values
        std::vector<real_t> values_;
    };

    void vec_copy(const Vector &src, Vector &dest);
    void vec_axpy(real_t a, const Vector &x, Vector &y);
}