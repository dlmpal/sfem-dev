#pragma once

namespace sfem::la
{
    /// @brief Linear algebra backend
    enum class Backend
    {
        native,
        petsc
    };
}