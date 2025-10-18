#include "fv_coefficient.hpp"

namespace sfem::fvm
{
    //=============================================================================
    FVCoefficient::FVCoefficient(std::shared_ptr<FVField> phi)
        : FieldCoefficient(phi), fv_phi_(phi)
    {
    }
    //=============================================================================
    std::shared_ptr<FVField> FVCoefficient::fv_field()
    {
        return fv_phi_;
    }
    //=============================================================================
    std::shared_ptr<const FVField> FVCoefficient::fv_field() const
    {
        return fv_phi_;
    }
    //=============================================================================
    FVCoefficient create_coeff(std::shared_ptr<const FVSpace> V,
                               const std::vector<std::string> &components)
    {
        auto phi = std::make_shared<FVField>(V, components);
        return FVCoefficient(phi);
    }
}