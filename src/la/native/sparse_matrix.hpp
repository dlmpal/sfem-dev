#pragma once

#include "../../graph/connectivity.hpp"
#include "../../parallel/index_map.hpp"
#include "../../base/config.hpp"
#include <memory>

// Forward declaration
namespace sfem::la
{
    class Vector;
}

namespace sfem::la
{
    /// @brief MPI-parallel sparse matrix in CSR format
    class SparseMatrix
    {
    public:
        /// @brief Create a SparseMatrix
        /// @param row_to_col Row-to-column connectivity
        /// @param row_index_map Row index map
        /// @param col_index_map Column index map
        /// @param block_size Block size (common for rows and columns)
        SparseMatrix(std::shared_ptr<const graph::Connectivity> row_to_col,
                     std::shared_ptr<const IndexMap> row_index_map,
                     std::shared_ptr<const IndexMap> col_index_map,
                     int block_size);

        /// @brief Get the row-to-column connectivity
        std::shared_ptr<const graph::Connectivity>
        connectivity() const;

        /// @brief Get the row and column index maps
        std::array<std::shared_ptr<const IndexMap>, 2>
        index_maps() const;

        /// @brief Get a reference to the matrix's values
        std::vector<real_t> &values();

        /// @brief Get a reference to the matrix's values (const version)
        const std::vector<real_t> &values() const;

        /// @brief Get the block size
        int block_size() const;

        /// @brief Set all matrix values to a uniform value
        void set_all(real_t value);

        /// @brief Set matrix values for a given set of (local) row and column indices
        /// @param row_idxs Local row indices
        /// @param col_idxs Local column indices
        /// @param values Values
        void set_values(std::span<const int> row_idxs,
                        std::span<const int> col_idxs,
                        std::span<const real_t> values);

        /// @brief Get the (local) column indices and a values for a given row
        /// @param row_idx Local row index
        std::pair<std::span<const int>, std::span<real_t>>
        row_data(int row_idx);

        /// @brief Get the (local) column indices and a values for a given row (const version)
        /// @param row_idx Local row index
        std::pair<std::span<const int>, std::span<const real_t>>
        row_data(int row_idx) const;

        /// @brief Assemble the matrix
        void assemble();

        /// @brief Fill an existing vector with the diagonal entries of the matrix
        void diagonal(Vector &diag) const;

        /// @brief Fill a component of an existing vector with the diagonal entries of the matrix
        /// for a given component
        void diagonal(Vector &diag, int src_comp, int dest_comp) const;

        /// @brief Scale the diagonal entries of the matrix
        void scale_diagonal(real_t a);

    private:
        /// @brief Row-to-column connectivity
        std::shared_ptr<const graph::Connectivity> row_to_col_;

        /// @brief Row index map
        std::shared_ptr<const IndexMap> row_im_;

        /// @brief Column index map
        std::shared_ptr<const IndexMap> col_im_;

        /// @brief Matrix values
        std::vector<real_t> values_;

        /// @brief Block size
        int bs_;
    };

    /// @brief Compute the Frobenius norm for a matrix
    real_t norm(const SparseMatrix &mat);

    /// @brief Sparse matrix-vector multiplication: y = Ax
    /// @note The ghost index values of x should be updated before calling
    void spmv(const SparseMatrix &A, const Vector &x, Vector &y);
}