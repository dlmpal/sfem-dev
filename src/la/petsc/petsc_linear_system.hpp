#pragma once

#ifdef SFEM_HAS_PETSC

#include "petsc_mat.hpp"
#include "petsc_vec.hpp"
#include "petsc_ksp.hpp"
#include "../native/linear_system.hpp"

namespace sfem::la::petsc
{
    class PetscLinearSystem : public LinearSystem
    {
    public:
        PetscLinearSystem(const IndexMap &index_map,
                          const graph::Connectivity &connectivity,
                          SolverType solver_type = SolverType::gmres,
                          SolverOptions solver_options = {},
                          int block_size = 1);

        PetscMat &A();
        const PetscMat &A() const;

        PetscVec &b();
        const PetscVec &b() const;

        PetscKSP &solver();
        const PetscKSP &solver() const;

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

    private:
        PetscMat A_;
        PetscVec b_;
        PetscKSP solver_;
        PetscVec diag_;
    };
}

#endif // SFEM_HAS_PETSC