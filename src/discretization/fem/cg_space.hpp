#pragma once

#include "fe_space.hpp"

namespace sfem::fem
{
    /// @brief Continuous Galerkin (CG) finite element space
    class CGSpace : public FESpace
    {
    public:
        CGSpace(std::shared_ptr<const mesh::Mesh> mesh, int order);
    };
}