#ifdef SFEM_HAS_PETSC

#include "petsc.hpp"
#include "petsc_vec.hpp"
#include "petsc_mat.hpp"
#include "petsc_ksp.hpp"
#include "../native/sparsity.hpp"
#include "../native/vector.hpp"
#include "../../base/logging.hpp"
#include <format>

namespace sfem::la::petsc
{
    //=============================================================================
    void initialize(int argc, char *argv[])
    {
        PetscInitialize(&argc, &argv, nullptr, nullptr);
    }
    //=============================================================================
    void finalize()
    {
        PetscFinalize();
    }
    //=============================================================================
    ISLocalToGlobalMapping create_is_from_im(const IndexMap &index_map, int block_size)
    {
        ISLocalToGlobalMapping is;
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,
                                     block_size,
                                     index_map.n_local(),
                                     index_map.local_idxs().data(),
                                     PETSC_COPY_VALUES,
                                     &is);
        return is;
    }
    //=============================================================================
    PetscVec create_vec(const IndexMap &index_map, int block_size)
    {
        Vec vec;
        if (block_size == 1)
        {
            VecCreateGhost(PETSC_COMM_WORLD,
                           index_map.n_owned(),
                           index_map.n_global(),
                           index_map.n_ghost(),
                           index_map.ghost_idxs().data(),
                           &vec);
        }
        else
        {
            VecCreateGhostBlock(PETSC_COMM_WORLD,
                                block_size,
                                index_map.n_owned() * block_size,
                                index_map.n_global() * block_size,
                                index_map.n_ghost(),
                                index_map.ghost_idxs().data(),
                                &vec);
        }

        // Set local-to-global mapping
        auto is = create_is_from_im(index_map, block_size);
        VecSetLocalToGlobalMapping(vec, is);

        return PetscVec(vec, true);
    }
    //=============================================================================
    PetscVec create_vec(const Vector &vec)
    {
        // Quick access
        const int block_size = vec.block_size();
        const auto index_map = vec.index_map();

        Vec vec_;
        if (vec.block_size() == 1)
        {
            VecCreateGhostWithArray(PETSC_COMM_WORLD,
                                    index_map->n_owned(),
                                    index_map->n_global(),
                                    index_map->n_ghost(),
                                    index_map->ghost_idxs().data(),
                                    vec.values().data(),
                                    &vec_);
        }
        else
        {
            VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,
                                         block_size,
                                         index_map->n_owned() * block_size,
                                         index_map->n_global() * block_size,
                                         index_map->n_ghost(),
                                         index_map->ghost_idxs().data(),
                                         vec.values().data(),
                                         &vec_);
        }

        // Set local-to-global mapping
        auto is = create_is_from_im(*index_map, block_size);
        VecSetLocalToGlobalMapping(vec_, is);

        return PetscVec(vec_, true);
    }
    //=============================================================================
    PetscMat create_mat(const graph::Connectivity &row_to_col,
                        const IndexMap &row_index_map,
                        const IndexMap &col_index_map,
                        int block_size)
    {
        auto [diag_nnz, off_diag_nnz] = compute_sparsity(row_to_col,
                                                         row_index_map,
                                                         col_index_map);
        Mat mat;
        if (block_size == 1)
        {
            MatCreateAIJ(PETSC_COMM_WORLD,
                         row_index_map.n_owned(),
                         col_index_map.n_owned(),
                         row_index_map.n_global(),
                         col_index_map.n_global(),
                         PETSC_DECIDE, ///< ignored
                         diag_nnz.data(),
                         PETSC_DECIDE, ///< ignored
                         off_diag_nnz.data(),
                         &mat);
        }
        else
        {
            MatCreateBAIJ(PETSC_COMM_WORLD,
                          block_size,
                          row_index_map.n_owned() * block_size,
                          col_index_map.n_owned() * block_size,
                          row_index_map.n_global() * block_size,
                          col_index_map.n_global() * block_size,
                          PETSC_DECIDE, ///< ignored
                          diag_nnz.data(),
                          PETSC_DECIDE, ///< ignored
                          off_diag_nnz.data(),
                          &mat);
        }

        ISLocalToGlobalMapping row_mapping;
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,
                                     block_size,
                                     row_index_map.n_local(),
                                     row_index_map.local_idxs().data(),
                                     PETSC_COPY_VALUES,
                                     &row_mapping);

        ISLocalToGlobalMapping col_mapping;
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,
                                     block_size,
                                     col_index_map.n_local(),
                                     col_index_map.local_idxs().data(),
                                     PETSC_COPY_VALUES,
                                     &col_mapping);

        MatSetLocalToGlobalMapping(mat, row_mapping, col_mapping);

        return PetscMat(mat, true);
    }
    //=============================================================================
    VecSet create_vecset(la::petsc::PetscVec &vec)
    {
        return [&vec](std::span<const int> idxs,
                      std::span<const real_t> values)
        {
            vec.set_values(idxs, values);
        };
    }
    //=============================================================================
    MatSet create_matset(la::petsc::PetscMat &mat)
    {
        return [&mat](std::span<const int> row_idxs,
                      std::span<const int> col_idxs,
                      std::span<const real_t> values)
        {
            mat.set_values(row_idxs, col_idxs, values);
        };
    }
    //=============================================================================
    void eliminate_rows_cols(std::span<const int> idxs,
                             std::span<const PetscReal> values,
                             PetscMat &A,
                             PetscVec &b,
                             PetscVec &x)
    {
        x.set_values(idxs, values, INSERT_VALUES);
        MatZeroRowsColumnsLocal(A.mat(),
                                static_cast<int>(idxs.size()),
                                idxs.data(), 1.0,
                                x.vec(), b.vec());
        x.assemble();
    }
    //=============================================================================
    void eliminate_rows_cols(std::span<const int> idxs, PetscMat &A)
    {
        MatZeroRowsColumnsLocal(A.mat(),
                                static_cast<int>(idxs.size()),
                                idxs.data(), 1.0,
                                nullptr, nullptr);
    }
    //=============================================================================
    int solve(const PetscMat &A,
              const PetscVec &b,
              PetscVec &x)
    {
        // Create and setup the linear solver
        la::petsc::PetscKSP solver;
        solver.set_from_options();
        solver.set_operator(A);

        // Solve and return iteration number;
        int n_iter = solver.solve(b, x);
        return n_iter;
    }
}

#endif // SFEM_USE_PETSC