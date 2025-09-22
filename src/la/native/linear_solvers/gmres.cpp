#include "gmres.hpp"
#include "../sparse_matrix.hpp"
#include "../../../base/error.hpp"

namespace sfem::la
{
    //=============================================================================
    GMRES::GMRES(real_t tol, int n_iter_max, bool verbose, int n_restart)
        : LinearSolver("GMRES", tol, n_iter_max, verbose),
          n_restart_(n_restart),
          v1(std::make_shared<IndexMap>(), 1),
          v2(std::make_shared<IndexMap>(), 1),
          x0(std::make_shared<IndexMap>(), 1),
          Q(1, 1), H(1, 1), e1(1, 1)
    {
    }
    //=============================================================================
    // Helper function: Set a dense matrix column
    // from the locally owned values of a vector
    static void set_col(DenseMatrix &Q,
                        const Vector &q, int idx)
    {
        SFEM_CHECK_SIZES(Q.n_rows(),
                         q.n_owned() * q.block_size());
        SFEM_CHECK_INDEX(idx, Q.n_cols());
        const int bs = q.block_size();
        for (int i = 0; i < q.n_owned(); i++)
        {
            for (int j = 0; j < bs; j++)
            {
                Q(i * bs + j, idx) = q(i, j);
            }
        }
    };
    //=============================================================================
    // Helper function: Set the locally owned values
    // of a vector from a column of a dense matrix
    static void get_col(const DenseMatrix &Q,
                        Vector &q, int idx)
    {
        SFEM_CHECK_SIZES(Q.n_rows(),
                         q.n_owned() * q.block_size());
        SFEM_CHECK_INDEX(idx, Q.n_cols());
        const int bs = q.block_size();
        for (int i = 0; i < q.n_owned(); i++)
        {
            for (int j = 0; j < bs; j++)
            {
                q(i, j) = Q(i * bs + j, idx);
            }
        }
    };
    //=============================================================================
    void GMRES::init(const SparseMatrix &A,
                     const Vector &b, Vector &x)
    {
        // Number of (local) equations
        const int m = b.n_owned() * b.block_size();

        // Degree of Krylov subspace
        const int n = n_restart_;

        // Allocate workspace objects
        v1 = Vector(b.index_map(), b.block_size());
        v2 = Vector(b.index_map(), b.block_size());
        x0 = Vector(b.index_map(), b.block_size());
        Q = DenseMatrix(m, n + 1);
        H = DenseMatrix(n + 1, n);
        e1 = DenseMatrix(n + 1, 1);

        //
        init_(0, A, b, x);
    }
    //=============================================================================
    void GMRES::init_(int iter, const SparseMatrix &A,
                      const Vector &b, Vector &x)
    {
        v1.set_all(0.0);
        v2.set_all(0.0);
        x0.set_all(0.0);
        Q.set_all(0.0);
        H.set_all(0.0);
        e1.set_all(0.0);

        //
        copy(x, x0);
        x0.update_ghosts();

        // Compute initial residual and its norm
        spmv(A, x0, v1);              // v1 = A*x0
        axpbypc(1, -1, 0, b, v1, v2); // v2 = b - v1 = b - A*x0
        residual_history_[iter] = norm(v2, NormType::l2);

        // Normalize v2 to obtain q1
        scale(1.0 / residual_history_[iter], v2);

        //
        set_col(Q, v2, 0);

        //
        e1(0, 0) = residual_history_[iter];

        //
        riter_ = 0;
    }
    //=============================================================================
    extern "C" void dgels_(const char *trans, const int *m, const int *n, const int *nrhs,
                           double *A, const int *lda, double *B, const int *ldb,
                           double *work, const int *lwork, int *info);
    DenseMatrix lstsq(const DenseMatrix &A,
                      const DenseMatrix &b)
    {
        SFEM_CHECK_SIZES(A.n_rows(), b.n_rows());

        const char trans = 'N';
        const int m = A.n_rows();
        const int n = A.n_cols();
        const int nrhs = 1;
        const int lda = A.n_rows();
        const int ldb = std::max(A.n_rows(), A.n_cols());
        auto A_trans = A.transpose();

        int info;
        int lwork = -1;
        double wkopt;

        std::vector<real_t> x_values(m);
        std::copy(b.values().cbegin(),
                  b.values().cend(),
                  x_values.begin());

        // Work array optimal size query
        dgels_(&trans, &m, &n, &nrhs,
               A_trans.values().data(), &lda,
               x_values.data(), &ldb,
               &wkopt, &lwork, &info);
        lwork = static_cast<int>(wkopt);
        if (info != 0)
        {
            SFEM_ERROR(std::format("dgels_ returned with info={}\n", info));
        }

        // Actual computation
        std::vector<real_t> work_arr(lwork);
        dgels_(&trans, &m, &n, &nrhs,
               A_trans.values().data(), &lda,
               x_values.data(), &ldb,
               work_arr.data(), &lwork, &info);
        x_values.resize(n);

        return DenseMatrix(n, 1, std::move(x_values));
    }
    //=============================================================================
    void GMRES::single_iteration(int iter, const SparseMatrix &A,
                                 const Vector &b, Vector &x)
    {
        // Number of (local) equations
        const int m = b.n_owned() * b.block_size();

        // Degree of Krylov subspace
        const int n = n_restart_;

        // Iteration since last restart
        const int k = riter_;

        // Perform a single Arnoldi iteration
        get_col(Q, v1, k); // v1 = q_(k-1)
        v1.update_ghosts();
        spmv(A, v1, v2); // v2 = A * v1 = A* q_(k-1)
        for (int j = 0; j < k + 1; j++)
        {
            get_col(Q, v1, j);      // v1 = q_j
            H(j, k) = dot(v1, v2);  // H_(j, k-1) = <v1, v2>
            axpy(-H(j, k), v1, v2); // v2 -= <v1, v2> v1
        }
        H(k + 1, k) = norm(v2, NormType::l2); // H(k, k-1) = <v2, v2>
        if (std::abs(H(k + 1, k)) >= std::numeric_limits<real_t>::epsilon() and k + 1 < n)
        {
            scale(1.0 / H(k + 1, k), v2); // v2 = v2 / ||v2||
        }
        set_col(Q, v2, k + 1);

        // Solve the least squares problem
        auto Hk = submatrix(H, 0, k + 2, 0, k + 1);
        auto e1k = submatrix(e1, 0, k + 2, 0, 1);
        auto yk = lstsq(Hk, e1k);

        // xk = x0 + Qk * yk
        const auto &x0_values = x0.values();
        auto &x_values = x.values();
        for (int i = 0; i < m; i++)
        {
            x_values[i] = x0_values[i];
            for (int j = 0; j < yk.n_rows(); j++)
            {
                x_values[i] += yk(j, 0) * Q(i, j);
            }
        }

        // Increment iteration (since last restart)
        riter_++;

        // Perform a restart, if required
        if (riter_ == n_restart_)
        {
            init_(iter, A, b, x);
        }
        else
        {
            // Update residual
            residual_history_[iter] = norm(Hk * yk - e1k);
        }
    }
}