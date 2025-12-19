#include "setval_utils.hpp"
#include <sfem/la/native/vector.hpp>
#include <sfem/la/native/sparse_matrix.hpp>

namespace sfem::la
{
    //=============================================================================
    VecSet create_vecset(Vector &vec)
    {
        return [&vec](std::span<const int> idxs,
                      std::span<const real_t> values)
        {
            vec.set_values(idxs, values);
        };
    }
    //=============================================================================
    MatSet create_matset(SparseMatrix &mat)
    {
        return [&mat](std::span<const int> row_idxs,
                      std::span<const int> col_idxs,
                      std::span<const real_t> values)
        {
            mat.set_values(row_idxs, col_idxs, values);
        };
    }
}
