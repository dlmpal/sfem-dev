#include "cg.hpp"
#include "../sparse_matrix.hpp"

namespace sfem::la
{
    //=============================================================================
    CG::CG(real_t tol, int n_iter_max, bool verbose)
        : LinearSolver("CG", tol, n_iter_max, verbose),
          Ap(std::make_shared<IndexMap>(), 1),
          p(std::make_shared<IndexMap>(), 1),
          r(std::make_shared<IndexMap>(), 1)
    {
    }
    //=============================================================================
    void CG::init(const SparseMatrix &A,
                  const Vector &b, Vector &x)
    {
        Ap = Vector(x.index_map(), x.block_size());
        spmv(A, x, Ap);

        r = Vector(x.index_map(), x.block_size());
        axpbypc(1, -1, 0, b, Ap, r);
        residual_history_[0] = norm(r, NormType::l2);

        p = Vector(x.index_map(), x.block_size());
        copy(r, p);
    }
    //=============================================================================
    void CG::single_iteration(int iter, const SparseMatrix &A,
                              [[maybe_unused]] const Vector &b, Vector &x)
    {
        // Compute Ap (intermediate product)
        p.update_ghosts();
        spmv(A, p, Ap);

        const real_t res_old = residual_history_[iter - 1];

        // Compute step size
        const real_t alpha = res_old * res_old / (dot(p, Ap));

        // Update solution vector x = x + alpha * p
        axpy(alpha, p, x);

        // Update residual vector: r = r - alpha * Ap
        axpy(-alpha, Ap, r);
        const real_t res_new = norm(r, NormType::l2);
        residual_history_[iter] = res_new;

        // Update search direction vector: p = r + beta * p
        const real_t beta = (res_new * res_new) / (res_old * res_old);
        axpbypc(1, beta, 0, r, p, p);
    }
}