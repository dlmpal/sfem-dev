#pragma once

#include <sfem/la/native/linear_solvers/linear_solver_factory.hpp>
#include <sfem/la/native/setval_utils.hpp>
#include <sfem/la/native/sparse_matrix.hpp>
#include <sfem/la/native/vector.hpp>

namespace sfem::la
{
    class LinearSystem
    {
    public:
        virtual ~LinearSystem() = default;

        virtual void reset() = 0;

        virtual MatSet lhs() = 0;
        virtual VecSet rhs() = 0;

        virtual void assemble() = 0;

        virtual void diagonal(la::Vector &diag) const = 0;
        virtual void scale_diagonal(real_t a) = 0;

        virtual void rhs_axpy(real_t a, const la::Vector &x) = 0;

        virtual void eliminate_dofs(std::span<const int> idxs,
                                    std::span<const real_t> values) = 0;

        virtual bool solve(Vector &x) = 0;

        virtual std::vector<real_t> residual_history() const = 0;
    };

    class NativeLinearSystem : public LinearSystem
    {
    public:
        NativeLinearSystem(std::shared_ptr<const IndexMap> index_map,
                           std::shared_ptr<const graph::Connectivity> connectivity,
                           SolverType solver_type = SolverType::gmres,
                           SolverOptions solver_options = {},
                           int block_size = 1);

        SparseMatrix &A();
        const SparseMatrix &A() const;

        Vector &b();
        const Vector &b() const;

        std::shared_ptr<LinearSolver> solver();
        std::shared_ptr<const LinearSolver> solver() const;

        void reset() override;

        MatSet lhs() override;
        VecSet rhs() override;

        void assemble() override;

        void diagonal(la::Vector &diag) const override;
        void scale_diagonal(real_t a) override;

        void rhs_axpy(real_t a, const la::Vector &x) override;

        void eliminate_dofs(std::span<const int> idxs,
                            std::span<const real_t> values) override;

        bool solve(Vector &x) override;

        std::vector<real_t> residual_history() const override;

    private:
        SparseMatrix A_;

        Vector b_;

        std::shared_ptr<LinearSolver> solver_;
    };
}