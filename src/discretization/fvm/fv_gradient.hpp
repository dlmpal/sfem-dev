#pragma once

#include "fv_bc.hpp"

namespace sfem::fvm
{
    enum GradientMethod
    {
        green_gauss = 0
    };

    ///
    void gradient(const FVFunction &phi,
                  const FVBC &bc,
                  FVFunction &grad,
                  GradientMethod method = GradientMethod::green_gauss);

    ///
    std::shared_ptr<FVFunction> gradient(const FVFunction &phi,
                                         const FVBC &bc,
                                         GradientMethod method = GradientMethod::green_gauss);

}