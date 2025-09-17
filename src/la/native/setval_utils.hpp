#pragma once

#include "../../base/config.hpp"
#include <functional>
#include <span>

// Forward declarations
namespace sfem::la
{
    class Vector;
    class SparseMatrix;
}

namespace sfem::la
{
    /// @brief Whether to set values in a buffer
    /// by inserting them (thereby overriding previously
    /// existing ones) or adding them
    enum class SetMode
    {
        add = 0,
        insert = 1
    };

    using VecSet = std::function<void(std::span<const int>,
                                      std::span<const real_t>)>;

    using MatSet = std::function<void(std::span<const int>,
                                      std::span<const int>,
                                      std::span<const real_t>)>;

    /// @brief Create a VecSet for a Vector
    VecSet create_vecset(Vector &vec);

    /// @brief Create a MatSet for a SparseMatrix
    MatSet create_matset(SparseMatrix &mat);
}