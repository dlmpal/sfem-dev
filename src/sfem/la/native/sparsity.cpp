#include "sparsity.hpp"
#include <sfem/parallel/mpi.hpp>

namespace sfem::la
{
    std::array<std::vector<int>, 2>
    compute_sparsity(const graph::Connectivity &row_to_col,
                     const IndexMap &row_index_map,
                     const IndexMap &col_index_map)
    {
        // Quick access
        int proc_rank = mpi::rank();

        // Count the number of non-zeros per row (for locally owned rows)
        std::vector<int> diag_nnz(row_index_map.n_owned(), 0);
        std::vector<int> off_diag_nnz(row_index_map.n_owned(), 0);
        for (int i = 0; i < row_to_col.n_primary(); i++)
        {
            int owner_i = row_index_map.get_owner(i);

            // Row "i" is locally owned
            if (owner_i == proc_rank)
            {
                // Loop over this row's columns
                for (int j : row_to_col.links(i))
                {
                    // Column "j" is locally owned
                    if (col_index_map.get_owner(j) == proc_rank)
                    {
                        diag_nnz[i]++;
                    }
                    else
                    {
                        off_diag_nnz[i]++;
                    }
                }
            }
        }

        return {diag_nnz, off_diag_nnz};
    }
}