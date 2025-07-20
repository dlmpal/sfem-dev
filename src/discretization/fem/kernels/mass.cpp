#include "mass.hpp"

namespace sfem::fem::kernels
{
    //=============================================================================
    MassND::MassND(int n_comp, std::shared_ptr<const Coefficient> coeff)
        : n_comp_(n_comp), coeff_(coeff)
    {
    }
    //=============================================================================
    la::DenseMatrix MassND::operator()(int cell_idx, const FEData &data)
    {
        const real_t coeff = (*coeff_)(cell_idx, 0);
        la::DenseMatrix M(data.N.n_rows() * n_comp_, data.N.n_rows() * n_comp_);
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            for (int j = 0; j < data.N.n_rows(); j++)
            {
                const real_t val = coeff * data.N(i, 0) * data.N(j, 0);
                for (int k = 0; k < n_comp_; k++)
                {
                    M(i * n_comp_ + k, j * n_comp_ + k) = val;
                }
            }
        }
        return M;
    }
}