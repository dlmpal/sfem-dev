#include "gmres.hpp"
#include "../sparse_matrix.hpp"
#include "../../../base/error.hpp"

namespace sfem::la
{
    //=============================================================================
    GMRES::GMRES(SolverOptions options, int n_restart)
        : LinearSolver("GMRES", options),
          n_restart_(n_restart),
          x0_(std::make_shared<IndexMap>(), 1),
          H_(n_restart_ + 1, n_restart_),
          e1_(n_restart_ + 1, 1)
    {
    }
    //=============================================================================
    void GMRES::init(const SparseMatrix &A,
                     const Vector &b, Vector &x)
    {
        // Allocate workspace objects
        x0_ = Vector(b.index_map(), b.block_size());
        Q_.clear();
        for (int i = 0; i < n_restart_ + 1; i++)
        {
            Q_.emplace_back(b.index_map(), b.block_size());
        }

        // Perform initial "restart"
        restart(0, A, b, x);
    }
    //=============================================================================
    void GMRES::restart(int iter, const SparseMatrix &A, const Vector &b, Vector &x)
    {
        // Save initial solution vector
        copy(x, x0_);
        x0_.update_ghosts();

        // Reset basis vectors
        for (auto &q : Q_)
        {
            q.set_all(0.0);
        }

        // Compute initial residual vector and its norm
        spmv(A, x0_, Q_[0]);
        axpy(-1, b, Q_[0]);
        residual_history_[iter] = norm(Q_[0], NormType::l2);

        // Normalize the initial residual vector to create the first basis vector
        // The negative sign is required since q0 was defined as Ax - b (instead of b - Ax)
        scale(-1.0 / residual_history_[iter], Q_[0]);

        // Reset Hessenberg matrix and e1 vector
        H_.set_all(0.0);
        e1_(0, 0) = residual_history_[iter];

        // Reset iterations since last restart
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
        const int k = riter_;

        // Perform a single Arnoldi iteration
        Q_[k].update_ghosts();
        spmv(A, Q_[k], Q_[k + 1]);
        for (int j = 0; j < k + 1; j++)
        {
            H_(j, k) = dot(Q_[j], Q_[k + 1]);
            axpy(-H_(j, k), Q_[j], Q_[k + 1]);
        }
        H_(k + 1, k) = norm(Q_[k + 1], NormType::l2);
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        if (std::abs(H_(k + 1, k)) >= eps and k + 1 < n_restart_)
        {
            scale(1.0 / H_(k + 1, k), Q_[k + 1]);
        }

        // Solve the least squares problem
        auto Hk = submatrix(H_, 0, k + 2, 0, k + 1);
        auto e1k = submatrix(e1_, 0, k + 2, 0, 1);
        // auto yk = lstsq(Hk, e1k);
        auto yk = Hk.invert().first * e1k;

        // Update solution vector
        copy(x0_, x);
        for (int j = 0; j < k + 1; j++)
        {
            axpy(yk(j, 0), Q_[j], x);
        }

        // Increment iteration (since last restart)
        riter_++;

        // Perform a restart, if required
        // Else, update residual history
        if (riter_ == n_restart_)
        {
            restart(iter, A, b, x);
        }
        else
        {
            residual_history_[iter] = norm(Hk * yk - e1k);
        }
    }
}