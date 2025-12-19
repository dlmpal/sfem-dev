#include "sparse_matrix.hpp"
#include <sfem/la/native/vector.hpp>
#include <sfem/parallel/mpi.hpp>
#include <sfem/base/error.hpp>
#include <cmath>
#include <numeric>

namespace sfem::la
{
    //=============================================================================
    SparseMatrix::SparseMatrix(std::shared_ptr<const graph::Connectivity> row_to_col,
                               std::shared_ptr<const IndexMap> row_index_map,
                               std::shared_ptr<const IndexMap> col_index_map,
                               int block_size)
        : row_to_col_(row_to_col),
          row_im_(row_index_map),
          col_im_(col_index_map),
          values_(row_to_col_->n_links() * block_size * block_size, 0.0),
          bs_(block_size)
    {
        SFEM_CHECK_SIZES(row_to_col_->n_primary(), row_im_->n_local());
        SFEM_CHECK_SIZES(row_to_col_->n_secondary(), col_im_->n_local());
    }
    //=============================================================================
    std::shared_ptr<const graph::Connectivity>
    SparseMatrix::connectivity() const
    {
        return row_to_col_;
    }
    //=============================================================================
    std::array<std::shared_ptr<const IndexMap>, 2>
    SparseMatrix::index_maps() const
    {
        return {row_im_, col_im_};
    }
    //=============================================================================
    std::vector<real_t> &SparseMatrix::values()
    {
        return values_;
    }
    //=============================================================================
    const std::vector<real_t> &SparseMatrix::values() const
    {
        return values_;
    }
    //=============================================================================
    int SparseMatrix::block_size() const
    {
        return bs_;
    }
    //=============================================================================
    void SparseMatrix::set_all(real_t value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }
    //=============================================================================
    void SparseMatrix::set_values(std::span<const int> row_idxs,
                                  std::span<const int> col_idxs,
                                  std::span<const real_t> values)
    {
        const int nr = static_cast<int>(row_idxs.size());
        const int nc = static_cast<int>(col_idxs.size());
        SFEM_CHECK_SIZES(nr * nc * bs_ * bs_, values.size());
        for (int i = 0; i < nr; i++)
        {
            const int r = row_idxs[i];
            const int ri = i * nc * bs_ * bs_;
            const int offset = row_to_col_->offset(r);
            for (int j = 0; j < nc; j++)
            {
                const int c = col_idxs[j];
                const int ci = j * bs_;
                const int rel_idx = row_to_col_->relative_index(r, c);
                const int start = (offset + rel_idx) * bs_ * bs_;
                for (int k1 = 0; k1 < bs_; k1++)
                {
                    for (int k2 = 0; k2 < bs_; k2++)
                    {
                        values_[start + k1 * bs_ + k2] += values[ri + ci + k1 * nc * bs_ + k2];
                    }
                }
            }
        }
    }
    //=============================================================================
    std::pair<std::span<const int>, std::span<real_t>>
    SparseMatrix::row_data(int row_idx)
    {
        return {row_to_col_->links(row_idx),
                {values_.begin() + row_to_col_->offset(row_idx) * bs_ * bs_,
                 values_.begin() + row_to_col_->offset(row_idx + 1) * bs_ * bs_}};
    }
    //=============================================================================
    std::pair<std::span<const int>, std::span<const real_t>>
    SparseMatrix::row_data(int row_idx) const
    {
        return {row_to_col_->links(row_idx),
                {values_.cbegin() + row_to_col_->offset(row_idx) * bs_ * bs_,
                 values_.cbegin() + row_to_col_->offset(row_idx + 1) * bs_ * bs_}};
    }
    //=============================================================================
    void SparseMatrix::assemble()
    {
        // Number of blocks for all ghost rows
        const int n_ghost_blocks = row_to_col_->offset(row_im_->n_local()) -
                                   row_to_col_->offset(row_im_->n_owned());

        // Save the ghost data in COO format
        std::vector<int> ghost_rows(n_ghost_blocks);
        std::vector<int> ghost_cols(n_ghost_blocks);
        std::vector<int> ghost_dest(n_ghost_blocks);
        std::vector<real_t> ghost_values(n_ghost_blocks * bs_ * bs_);
        int displ = 0;
        for (int i = row_im_->n_owned(); i < row_im_->n_local(); i++)
        {
            std::fill(ghost_rows.begin() + displ,
                      ghost_rows.begin() + displ + row_to_col_->n_links(i),
                      i);

            std::fill(ghost_dest.begin() + displ,
                      ghost_dest.begin() + displ + row_to_col_->n_links(i),
                      row_im_->get_owner(i));

            auto [row_cols, row_values] = row_data(i);

            std::copy(row_cols.cbegin(),
                      row_cols.cend(),
                      ghost_cols.begin() + displ);

            std::copy(row_values.cbegin(),
                      row_values.cend(),
                      ghost_values.begin() + displ * bs_ * bs_);

            // Set all values of ghost rows to zero
            std::fill(row_values.begin(), row_values.end(), 0.0);

            displ += row_to_col_->n_links(i);
        }

        // Convert ghost rows and columns to global indexing
        ghost_rows = row_im_->local_to_global(ghost_rows);
        ghost_cols = col_im_->local_to_global(ghost_cols);

        // Send ghost data to destination (owner) processes
        auto recv_rows = std::get<0>(mpi::send_to_dest<int>(ghost_rows, ghost_dest));
        auto recv_cols = std::get<0>(mpi::send_to_dest<int>(ghost_cols, ghost_dest));
        auto recv_values = std::get<0>(mpi::send_to_dest<real_t>(ghost_values, ghost_dest, bs_ * bs_));

        // Convert received rows and columns to local indexing
        recv_rows = row_im_->global_to_local(recv_rows);
        recv_cols = col_im_->global_to_local(recv_cols);

        // Add the contributions
        for (std::size_t i = 0; i < recv_rows.size(); i++)
        {
            // Compute the start of the current received block
            const int offset = row_to_col_->offset(recv_rows[i]);
            const int rel_idx = row_to_col_->relative_index(recv_rows[i], recv_cols[i]);
            const int start = (offset + rel_idx) * bs_ * bs_;

            for (int k1 = 0; k1 < bs_; k1++)
            {
                for (int k2 = 0; k2 < bs_; k2++)
                {
                    values_[start + k1 * bs_ + k2] += recv_values[i * bs_ * bs_ + k1 * bs_ + k2];
                }
            }
        }
    }
    //=============================================================================
    void SparseMatrix::diagonal(Vector &diag) const
    {
        SFEM_CHECK_SIZES(row_im_->n_owned(), diag.n_owned());
        SFEM_CHECK_SIZES(bs_, diag.block_size());
        for (int r = 0; r < row_im_->n_owned(); r++)
        {
            const int offset = row_to_col_->offset(r);
            const int rel_idx = row_to_col_->relative_index(r, r);
            const int start = (offset + rel_idx) * bs_ * bs_;
            for (int k = 0; k < bs_; k++)
            {
                diag(r, k) = values_[start + k * bs_ + k];
            }
        }
    }
    //=============================================================================
    void SparseMatrix::diagonal(Vector &diag, int src_comp, int dest_comp) const
    {
        SFEM_CHECK_SIZES(row_im_->n_owned(), diag.n_owned());
        SFEM_CHECK_INDEX(src_comp, bs_);
        for (int r = 0; r < row_im_->n_owned(); r++)
        {
            const int offset = row_to_col_->offset(r);
            const int rel_idx = row_to_col_->relative_index(r, r);
            const int start = (offset + rel_idx) * bs_ * bs_;
            diag(r, dest_comp) = values_[start + src_comp * bs_ + src_comp];
        }
    }
    //=============================================================================
    void SparseMatrix::scale_diagonal(real_t a)
    {
        for (int r = 0; r < row_im_->n_owned(); r++)
        {
            const int offset = row_to_col_->offset(r);
            const int rel_idx = row_to_col_->relative_index(r, r);
            const int start = (offset + rel_idx) * bs_ * bs_;
            for (int k = 0; k < bs_; k++)
            {
                values_[start + k * bs_ + k] *= a;
            }
        }
    }
    //=============================================================================
    real_t norm(const SparseMatrix &A)
    {
        real_t norm = std::accumulate(A.values().cbegin(),
                                      A.values().cend(), 0.0,
                                      [](real_t acc, real_t v)
                                      { return acc + v * v; });
        return std::sqrt(mpi::reduce(norm, mpi::ReduceOperation::sum));
    }
    //=============================================================================
    void spmv(const SparseMatrix &A,
              const Vector &x,
              Vector &y)
    {
        const int bs = A.block_size();
        const auto row_im = A.index_maps()[0];
        const auto col_im = A.index_maps()[1];

        SFEM_CHECK_SIZES(row_im->n_owned(), y.index_map()->n_owned());
        SFEM_CHECK_SIZES(col_im->n_owned(), x.index_map()->n_owned());
        SFEM_CHECK_SIZES(bs, x.block_size());
        SFEM_CHECK_SIZES(bs, y.block_size());

        y.set_all(0.0);
        for (int r = 0; r < row_im->n_owned(); r++)
        {
            const auto [cols, values] = A.row_data(r);
            for (std::size_t c = 0; c < cols.size(); c++)
            {
                for (int k1 = 0; k1 < bs; k1++)
                {
                    for (int k2 = 0; k2 < bs; k2++)
                    {
                        y(r, k1) += x(cols[c], k2) * values[c * bs * bs + k1 * bs + k2];
                    }
                }
            }
        }
    }
}