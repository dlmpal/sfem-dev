#include "diffusion.hpp"

namespace sfem::fem::kernels
{
    //=============================================================================
    Diffusion2D::Diffusion2D(std::shared_ptr<const Coefficient> coeff)
        : coeff_(coeff)
    {
    }
    //=============================================================================
    la::DenseMatrix Diffusion2D::operator()(int cell_idx, const fem::FEData &data)
    {
        const real_t coeff = (*coeff_)(cell_idx, 0);
        la::DenseMatrix K(data.N.n_rows(), data.N.n_rows());
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            for (int j = 0; j < data.N.n_rows(); j++)
            {
                K(i, j) = coeff * (data.dNdX(i, 0) * data.dNdX(j, 0) +
                                   data.dNdX(i, 1) * data.dNdX(j, 1));
            }
        }
        return K;
    }
    //=============================================================================
    Diffusion3D::Diffusion3D(std::shared_ptr<const Coefficient> coeff)
        : coeff_(coeff)
    {
    }
    //=============================================================================
    la::DenseMatrix Diffusion3D::operator()(int cell_idx, const fem::FEData &data)
    {
        const real_t coeff = (*coeff_)(cell_idx, 0);
        la::DenseMatrix K(data.N.n_rows(), data.N.n_rows());
        for (int i = 0; i < data.N.n_rows(); i++)
        {
            for (int j = 0; j < data.N.n_rows(); j++)
            {
                K(i, j) = coeff * (data.dNdX(i, 0) * data.dNdX(j, 0) +
                                   data.dNdX(i, 1) * data.dNdX(j, 1) +
                                   data.dNdX(i, 2) * data.dNdX(j, 2));
            }
        }
        return K;
    }
}