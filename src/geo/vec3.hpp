#pragma once

#include "../base/config.hpp"
#include <array>

namespace sfem::geo
{
    /// @brief 3D vector
    class Vec3
    {
    public:
        /// @brief Create a vector equal to (x, y, z)
        Vec3(real_t x = 0.0, real_t y = 0.0, real_t z = 0.0);

        /// @brief Create a vector from point-1 to point-2
        Vec3(const std::array<real_t, 3> &point1,
             const std::array<real_t, 3> &point2);

        /// @brief Get the vector's x component
        real_t x() const;

        /// @brief Get the vector's y component
        real_t y() const;

        /// @brief Get the vector's z component
        real_t z() const;

        /// @brief Get the vector's component in a given direction
        real_t operator()(int dir) const;

        /// @brief Add another vector
        Vec3 &operator+=(const Vec3 &rhs);

        /// @brief Add a scalar constant to the vector
        Vec3 &operator+=(real_t rhs);

        /// @brief Subtract another vector
        Vec3 &operator-=(const Vec3 &rhs);

        /// @brief Subtract a scalar constant from the vector
        Vec3 &operator-=(real_t rhs);

        /// @brief Scale the vector by a scalar constant
        Vec3 &operator*=(real_t rhs);

        /// @brief Divide the vector by a scalar constant
        Vec3 &operator/=(real_t rhs);

        /// @brief Get the vector's magnitude
        real_t mag() const;

        /// @brief Return the normalized vector
        Vec3 normalize() const;

        /// @brief Returns a vector normal to the original, in the 2D sense.
        /// Thus for a vector v = (x, y, z),
        /// the normal is defined as n = (y, -x, z)
        Vec3 normal2D() const;

    private:
        /// @brief Vector xyz components
        std::array<real_t, 3> x_;
    };

    /// @brief Cross product of two 3D vectors
    Vec3 cross(const Vec3 &lhs, const Vec3 &rhs);

    /// @brief Inner product of two 3D vectors
    real_t inner(const Vec3 &lhs, const Vec3 &rhs);

    /// @brief Scale a vector by a scalar constant
    Vec3 operator*(const Vec3 &lhs, real_t rhs);

    /// @brief Scale a vector by a scalar constant
    Vec3 operator*(real_t lhs, const Vec3 &rhs);

    /// @brief Divide a vector by a scalar constant
    Vec3 operator/(const Vec3 &lhs, real_t rhs);

    /// @brief Add a scalar constant to a vector
    Vec3 operator+(const Vec3 &lhs, real_t rhs);

    /// @brief Add a scalar constant to a vector
    Vec3 operator+(real_t lhs, const Vec3 &rhs);

    /// @brief Subtract a scalar constant from a vector
    Vec3 operator-(const Vec3 &lhs, real_t rhs);

    /// @brief Subtract a scalar constant from a vecto
    Vec3 operator-(real_t lhs, const Vec3 &rhs);

    /// @brief Add two vectors
    Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs);

    /// @brief Subtract a vector from another
    Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs);
}