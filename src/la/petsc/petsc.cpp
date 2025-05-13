#ifdef SFEM_HAS_PETSC

#include "petsc.hpp"
#include "petsc_vec.hpp"
#include "petsc_mat.hpp"
#include "petsc_ksp.hpp"
#include "../native/sparsity.hpp"
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
    void error(int error_code, std::source_location location)
    {
        const char *error_msg;
        PetscErrorMessage(error_code, &error_msg, nullptr);
        log_msg(std::format("Call to PETSc function returned with error:\n\t{}\n", error_msg),
                LogLevel::error, location);
    }
    //=============================================================================
    PetscVec create_vec(const IndexMap &index_map, int block_size)
    {
        Vec vec;
        int err_code;
        if (block_size == 1)
        {
            err_code = VecCreateGhost(PETSC_COMM_WORLD,
                                      index_map.n_owned(),
                                      index_map.n_global(),
                                      index_map.n_ghost(),
                                      index_map.ghost_idxs().data(),
                                      &vec);
        }
        else
        {
            err_code = VecCreateGhostBlock(PETSC_COMM_WORLD,
                                           block_size,
                                           index_map.n_owned() * block_size,
                                           index_map.n_global() * block_size,
                                           index_map.n_ghost(),
                                           index_map.ghost_idxs().data(),
                                           &vec);
        }
        SFEM_CHECK_PETSC_ERROR(err_code);

        ISLocalToGlobalMapping mapping;
        err_code = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,
                                                block_size,
                                                index_map.n_local(),
                                                index_map.local_idxs().data(),
                                                PETSC_COPY_VALUES,
                                                &mapping);
        SFEM_CHECK_PETSC_ERROR(err_code);

        err_code = VecSetLocalToGlobalMapping(vec, mapping);
        SFEM_CHECK_PETSC_ERROR(err_code);

        return PetscVec(vec, true);
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
        int err_code;
        if (block_size == 1)
        {
            err_code = MatCreateAIJ(PETSC_COMM_WORLD,
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
            err_code = MatCreateBAIJ(PETSC_COMM_WORLD,
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
        SFEM_CHECK_PETSC_ERROR(err_code);

        ISLocalToGlobalMapping row_mapping;
        err_code = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,
                                                block_size,
                                                row_index_map.n_local(),
                                                row_index_map.local_idxs().data(),
                                                PETSC_COPY_VALUES,
                                                &row_mapping);
        SFEM_CHECK_PETSC_ERROR(err_code);

        ISLocalToGlobalMapping col_mapping;
        err_code = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,
                                                block_size,
                                                col_index_map.n_local(),
                                                col_index_map.local_idxs().data(),
                                                PETSC_COPY_VALUES,
                                                &col_mapping);
        SFEM_CHECK_PETSC_ERROR(err_code);

        err_code = MatSetLocalToGlobalMapping(mat, row_mapping, col_mapping);
        SFEM_CHECK_PETSC_ERROR(err_code);

        return PetscMat(mat, true);
    }
    //=============================================================================
    void eliminate_rows_cols(std::span<const int> idxs,
                             std::span<const PetscReal> values,
                             PetscMat &A,
                             PetscVec &b,
                             PetscVec &x)
    {
        x.set_values(idxs, values, true);
        int error_code = MatZeroRowsColumnsLocal(A.mat(),
                                                 static_cast<int>(idxs.size()),
                                                 idxs.data(),
                                                 1.0, x.vec(), b.vec());
        SFEM_CHECK_PETSC_ERROR(error_code);
        x.assemble();
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