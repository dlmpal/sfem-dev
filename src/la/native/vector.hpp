#pragma once

#include "setval_utils.hpp"
#include "../../parallel/index_map.hpp"
#include <memory>

namespace sfem::la
{
    /// @brief Distributed vector the values of which are partitioned
    /// according to a given index map. The vector stores values for both
    /// locally owned and ghost indices.
    ///
    /// If vector values are inserted incrementally (e.g. in finite-element assembly)
    /// the user should first call assemble(), before trying to access vector values.
    /// Calling assemble() will zero-out ghost index values after sending them to their
    /// respective owner processes. If ghost index values are also to be accessed,
    /// the user should also call update_ghosts() immediately after assemble().
    ///
    /// Vector operations (e.g adding, scaling and copying vector values) are performed locally,
    /// i.e. ONLY for owned indices. As was the case for vectors the values of which where incrementally
    /// inserted, update_ghosts() should be called if ghost index values are to be accessed for vectors
    /// which are the result of vector operations.
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
        /// @note Includes ghost index values
        std::vector<real_t> &values();

        /// @brief Get the vector's values (const version)
        /// @note Includes ghost index values
        const std::vector<real_t> &values() const;

        /// @brief Get the number of owned indices
        int n_owned() const;

        /// @brief Get the number of ghost indices
        int n_ghost() const;

        /// @brief Get the number of local indices
        /// @note Sum of owned and ghost indices
        int n_local() const;

        /// @brief Get the global number of indices
        int n_global() const;

        /// @brief Get the value for a given (local) index and component
        real_t &operator()(int idx, int comp = 0);

        /// @brief Get the value for a given (local) index and component  (const version)
        real_t operator()(int idx, int comp = 0) const;

        /// @brief Set all vector values to a uniform value
        /// @note Includes ghost values
        void set_all(real_t value);

        /// @brief Set vector values for a given set of (local) indices
        /// @param idxs Local indices
        /// @param values Values
        /// @param mode Whether to insert or add the values
        void set_values(std::span<const int> idxs,
                        std::span<const real_t> values,
                        SetMode mode = SetMode::add);

        /// @brief Assemble the vector
        void assemble();

        /// @brief Update the values of ghost indices
        void update_ghosts();

    protected:
        /// @brief Index map
        std::shared_ptr<const IndexMap> im_;

        /// @brief Block size
        int bs_;

        /// @brief Vector values
        std::vector<real_t> values_;
    };

    // The following vector operations are performed locally,
    // i.e. for locally owned indices ONLY

    /// @brief Copy the values of one vector to another
    /// @note Values of ghost indices are NOT copied
    void copy(const Vector &src, Vector &dest);

    /// @brief Perform the operation x = a * x
    void scale(real_t a, Vector &x);

    /// @brief Perform the operation: y = y + a * x
    void axpy(real_t a, const Vector &x, Vector &y);

    /// @brief Perform the operation: z = a * x + b * y + c
    void axpbypc(real_t a, real_t b, real_t c,
                 const Vector &x, const Vector &y, Vector &z);

    /// @brief Compute the dot product for a pair of vectors
    real_t dot(const Vector &x, const Vector &y);

    enum class NormType
    {
        l1,
        l2,
        linf
    };

    /// @brief Compute a norm of a given type for a vector
    real_t norm(const Vector &x, NormType norm_type);
}