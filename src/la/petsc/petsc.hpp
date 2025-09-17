#pragma once

#ifdef SFEM_HAS_PETSC

#include "../native/setval_utils.hpp"
#include "../../graph/connectivity.hpp"
#include "../../parallel/index_map.hpp"
#include <petsc.h>
#include <source_location>

#if defined(SFEM_USE_SINGLE_PRECISION) and !defined(PETSC_USE_REAL_SINGLE)
#error "Mismatch between SFEM real type and PETSc real type"
#endif
#if defined(SFEM_USE_DOUBLE_PRECISION) and !defined(PETSC_USE_REAL_DOUBLE)
#error "Mismatch between SFEM real type and PETSc real type"
#endif
#ifdef PETSC_USE_COMPLEX
#error "SFEM cannot utilize complex arithmetic in PETSc"
#endif
#ifdef PETSC_USE_64_BIT_INDICES
#error "SFEM does not currently support 64-bit integers"
#endif

// Forward declarations
namespace sfem::la
{
    class Vector;
}

namespace sfem::la::petsc
{
    // Convenience functions for initializing and finalizing PETSc
    void initialize(int argc, char *argv[]);
    void finalize();

    /// @brief Create an ISLocalToGlobalMapping from an IndexMap and for a given block size
    ISLocalToGlobalMapping create_is_from_im(const IndexMap &index_map, int block_size);

    // Forward declarations
    // These classes are thin wrappers around PETSc's Vec, Mat
    // and KSP objects, with the main purpose of automatic memory
    // management
    class PetscVec;
    class PetscMat;
    class PetscKSP;

    /// @brief Create a PETSc Vec for a specified index map and block size
    PetscVec create_vec(const IndexMap &index_map, int block_size);

    /// @brief Create a PETSc Vec from a native Vector
    /// @note The underlying memory is owned by the Vector
    PetscVec create_vec(const Vector &vec);

    /// @brief Create a PETSc Mat for specified row-to-column connectivity,
    /// row and column index maps and block size
    PetscMat create_mat(const graph::Connectivity &row_to_col,
                        const IndexMap &row_index_map,
                        const IndexMap &col_index_map,
                        int block_size);

    /// @brief Create a VecSet for a PETSc vec
    VecSet create_vecset(la::petsc::PetscVec &vec);

    /// @brief Create a MatSet for a PETSc mat
    MatSet create_matset(la::petsc::PetscMat &mat);

    /// @brief For a linear system of the form Ax=b,
    /// remove rows and columns of the LHS (A) corresponding,
    /// for example, to essential boundary conditions and add
    /// their contributios to the RHS (b)
    /// @note Does not affect the sparsity of A
    /// @note Also calls assemble() on x
    /// @param idxs Indices of the fixed DoF
    /// @param values Values of the fixed DoF
    /// @param A Left-hand-side (LHS) matrix
    /// @param b Right-hand-side (RHS) vector
    /// @param x Solution vector
    void eliminate_rows_cols(std::span<const int> idxs,
                             std::span<const PetscReal> values,
                             PetscMat &A,
                             PetscVec &b,
                             PetscVec &x);

    /// @brief For a linear system of the form Ax=b,
    /// remove rows and columns of the LHS (A) corresponding,
    /// for example, to essential boundary conditions.
    /// @note Does not affect the sparsity of A
    /// @note Also calls assemble() on x
    /// @param idxs Indices of the fixed DoF
    /// @param A Left-hand-side (LHS) matrix
    void eliminate_rows_cols(std::span<const int> idxs, PetscMat &A);

    /// @brief Solve the linear system Ax=b using PETSc's KSP solvers
    /// @param A Left-hand-side (LHS) matrix
    /// @param b Right-hand-side (RHS) vector
    /// @param x Solution vector
    /// @note This creates a PetscKSP at each call
    /// @return Number of iterations performed
    int solve(const PetscMat &A,
              const PetscVec &b,
              PetscVec &x);
}

#endif // SFEM_HAS_PETSC