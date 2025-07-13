#pragma once

#include "../fe_space.hpp"
#include "../../../la/native/vector.hpp"

namespace sfem::fem
{
    using FECellKernel = std::function<la::DenseMatrix(int, const FEData &)>;
    using FEFacetKernel = std::function<la::DenseMatrix(int, const FEData &, const geo::Vec3 &)>;
    using VecSet = std::function<void(std::span<const int>, std::span<const real_t>)>;
    using MatSet = std::function<void(std::span<const int>, std::span<const int>, std::span<const real_t>)>;

    void assemble_matrix_cells(const FESpace &phi,
                               const std::string &region,
                               FECellKernel kernel,
                               MatSet mat);

    void assemble_matrix_facets(const FESpace &phi,
                                const std::string &region,
                                FEFacetKernel kernel,
                                MatSet mat);

    void assemble_vec_cells(const FESpace &phi,
                            const std::string &region,
                            FECellKernel kernel,
                            VecSet vec);

    void assemble_vec_facets(const FESpace &phi,
                             const std::string &region,
                             FEFacetKernel kernel,
                             VecSet vec);

    la::DenseMatrix integrate_cells(const FESpace &phi,
                                    const std::string &region,
                                    FECellKernel kernel);

    la::DenseMatrix integrate_facets(const FESpace &phi,
                                     const std::string &region,
                                     FEFacetKernel kernel);
}