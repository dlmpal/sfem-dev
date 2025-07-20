#pragma once

#include "../fe_function.hpp"
#include "../../../la/utils.hpp"

namespace sfem::fem
{
    using FECellKernel = std::function<la::DenseMatrix(int, const FEData &)>;
    using FEFacetKernel = std::function<la::DenseMatrix(int, const FEData &, const geo::Vec3 &)>;

    void assemble_matrix_cells(const FEFunction &phi, const std::string &region,
                               FECellKernel kernel, la::MatSet mat);

    void assemble_matrix_facets(const FEFunction &phi, const std::string &region,
                                FEFacetKernel kernel, la::MatSet mat);

    void assemble_vec_cells(const FEFunction &phi, const std::string &region,
                            FECellKernel kernel, la::VecSet vec);

    void assemble_vec_facets(const FEFunction &phi, const std::string &region,
                             FEFacetKernel kernel, la::VecSet vec);

    la::DenseMatrix integrate_cells(const FEFunction &phi, const std::string &region,
                                    FECellKernel kernel);

    la::DenseMatrix integrate_facets(const FEFunction &phi, const std::string &region,
                                     FEFacetKernel kernel);
}