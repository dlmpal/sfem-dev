#pragma once

#include "fv_field.hpp"
#include "../coefficient.hpp"

namespace sfem::fvm
{
    class FVCoefficient : public FieldCoefficient
    {
    public:
        FVCoefficient(std::shared_ptr<FVField> phi);

        std::shared_ptr<FVField> fv_field();
        std::shared_ptr<const FVField> fv_field() const;

    private:
        std::shared_ptr<FVField> fv_phi_;
    };

    FVCoefficient create_coeff(std::shared_ptr<const FVSpace> V,
                               const std::vector<std::string> &components);
}