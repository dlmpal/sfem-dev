#include "strain.hpp"

namespace sfem::fem::solid_mechanics
{
    //=============================================================================
    Strain::Strain(FEField U)
        : U_(U)
    {
        if (U_.n_comp() < 2 or U.n_comp() > 3)
        {
            SFEM_ERROR(std::format("Invalid field (n_comp={})\n", U.n_comp()));
        }
    }
    //=============================================================================
    FEField &Strain::U()
    {
        return U_;
    }
    //=============================================================================
    const FEField &Strain::U() const
    {
        return U_;
    }
    //=============================================================================
    int Strain::n_strain() const
    {
        if (U_.n_comp() == 2)
        {
            return 3;
        }
        else
        {
            return 6;
        }
    }
    //=============================================================================
    SmallStrain::SmallStrain(FEField U)
        : Strain(U)
    {
    }
    //=============================================================================
    void SmallStrain::B_geo(const FEData &, la::DenseMatrix &) const
    {
    }
    //=============================================================================
    void SmallStrain::B_mat(const FEData &data, la::DenseMatrix &B) const
    {
        for (int i = 0; i < data.n_nodes; i++)
        {
            if (data.pdim == 2)
            {
                B(0, i * 2 + 0) = data.dNdX(i, 0); ///< exx
                B(1, i * 2 + 1) = data.dNdX(i, 1); ///< eyy
                B(2, i * 2 + 0) = data.dNdX(i, 1); ///< exy
                B(2, i * 2 + 1) = data.dNdX(i, 0); ///< exy
            }
            else
            {
                B(0, i * 3 + 0) = data.dNdX(i, 0); ///< exx
                B(1, i * 3 + 1) = data.dNdX(i, 1); ///< eyy
                B(2, i * 3 + 2) = data.dNdX(i, 2); ///< ezz

                B(3, i * 3 + 0) = data.dNdX(i, 1); ///< exy
                B(3, i * 3 + 1) = data.dNdX(i, 0); ///< exy

                B(4, i * 3 + 1) = data.dNdX(i, 2); ///< eyz
                B(4, i * 3 + 2) = data.dNdX(i, 1); ///< eyz

                B(5, i * 3 + 0) = data.dNdX(i, 2); ///< ezx
                B(5, i * 3 + 2) = data.dNdX(i, 0); ///< ezx
            }
        }
    }
    //=============================================================================
    void SmallStrain::operator()(const FEData &data, la::DenseMatrix &e) const
    {
        const int n_dof = data.n_nodes * data.pdim;
        la::DenseMatrix u(n_dof, 1);
        U_.cell_values(data.elem_idx, u.values());
        la::DenseMatrix B(n_strain(), n_dof);
        B_mat(data, B);
        e = B * u;
    }
}