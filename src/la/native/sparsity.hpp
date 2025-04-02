#pragma once

#include "../../graph/connectivity.hpp"
#include "../../parallel/index_map.hpp"

namespace sfem::la
{
    std::array<std::vector<int>, 2>
    compute_sparsity(const graph::Connectivity &row_to_col,
                     const IndexMap &row_index_map,
                     const IndexMap &col_index_map);
}