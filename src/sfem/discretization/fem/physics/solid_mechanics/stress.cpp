#include "stress.hpp"

namespace sfem::fem::solid_mechanics
{
    //=============================================================================
    Stress::Stress(FEField U, Strain &strain, ElasticityConstitutive &constitutive)
        : U_(U),
          strain_(strain),
          constitutive_(constitutive)
    {
    }
    //=============================================================================
    FEField &Stress::U()
    {
        return U_;
    }
    //=============================================================================
    const FEField &Stress::U() const
    {
        return U_;
    }
    //=============================================================================
    Strain &Stress::strain()
    {
        return strain_;
    }
    //=============================================================================
    const Strain &Stress::strain() const
    {
        return strain_;
    }
    //=============================================================================
    ElasticityConstitutive &Stress::constitutive()
    {
        return constitutive_;
    }
    //=============================================================================
    const ElasticityConstitutive &Stress::constitutive() const
    {
        return constitutive_;
    }
    //=============================================================================
    void Stress::operator()(const FEData &data, la::DenseMatrix &s) const
    {
        la::DenseMatrix e(strain_.n_strain(), 1);
        strain_(data, e);
        constitutive_.stress(data, e, s);
    }
}