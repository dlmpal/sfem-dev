#include "vec3.hpp"
#include <cmath>

namespace sfem::geo
{
    //=============================================================================
    Vec3::Vec3(real_t x, real_t y, real_t z)
        : x_{x, y, z}
    {
    }
    //=============================================================================
    Vec3::Vec3(const std::array<real_t, 3> &point1,
               const std::array<real_t, 3> &point2)
        : Vec3(point2[0] - point1[0],
               point2[1] - point1[1],
               point2[2] - point1[2])
    {
    }
    //=============================================================================
    real_t Vec3::x() const
    {
        return x_[0];
    }
    //=============================================================================
    real_t Vec3::y() const
    {
        return x_[1];
    }
    //=============================================================================
    real_t Vec3::z() const
    {
        return x_[2];
    }
    //=============================================================================
    real_t &Vec3::operator()(int dir)
    {
        if (dir < 0 or dir > 2)
        {
            /// @todo error
        }

        return x_[dir];
    }
    //=============================================================================
    real_t Vec3::operator()(int dir) const
    {
        if (dir < 0 or dir > 2)
        {
            /// @todo error
        }

        return x_[dir];
    }
    //=============================================================================
    Vec3 &Vec3::operator+=(const Vec3 &rhs)
    {
        x_[0] += rhs.x_[0];
        x_[1] += rhs.x_[1];
        x_[2] += rhs.x_[2];
        return *this;
    }
    //=============================================================================
    Vec3 &Vec3::operator+=(real_t rhs)
    {
        x_[0] += rhs;
        x_[1] += rhs;
        x_[2] += rhs;
        return *this;
    }
    //=============================================================================
    Vec3 &Vec3::operator-=(const Vec3 &rhs)
    {
        x_[0] -= rhs.x_[0];
        x_[1] -= rhs.x_[1];
        x_[2] -= rhs.x_[2];
        return *this;
    }
    //=============================================================================
    Vec3 &Vec3::operator-=(real_t rhs)
    {
        x_[0] -= rhs;
        x_[1] -= rhs;
        x_[2] -= rhs;
        return *this;
    }
    //=============================================================================
    Vec3 &Vec3::operator*=(real_t rhs)
    {
        x_[0] *= rhs;
        x_[1] *= rhs;
        x_[2] *= rhs;
        return *this;
    }
    //=============================================================================
    Vec3 &Vec3::operator/=(real_t rhs)
    {
        rhs = 1 / rhs;
        x_[0] *= rhs;
        x_[1] *= rhs;
        x_[2] *= rhs;
        return *this;
    }
    //=============================================================================
    real_t Vec3::mag() const
    {
        return std::sqrt(x_[0] * x_[0] + x_[1] * x_[1] + x_[2] * x_[2]);
    }
    //=============================================================================
    Vec3 Vec3::normalize() const
    {
        real_t mag_ = mag();
        return Vec3(x_[0] / mag_, x_[1] / mag_, x_[2] / mag_);
    }
    //=============================================================================
    Vec3 Vec3::normal2D() const
    {
        return Vec3(x_[1], -x_[0], x_[2]);
    }
    //=============================================================================
    Vec3 cross(const Vec3 &lhs, const Vec3 &rhs)
    {
        return Vec3(lhs.y() * rhs.z() - lhs.z() * rhs.y(),
                    lhs.z() * rhs.x() - lhs.x() * rhs.z(),
                    lhs.x() * rhs.y() - lhs.y() * rhs.x());
    }
    //=============================================================================
    real_t inner(const Vec3 &lhs, const Vec3 &rhs)
    {
        return lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z();
    }
    //=============================================================================
    Vec3 operator*(const Vec3 &lhs, real_t rhs)
    {
        return Vec3(lhs.x() * rhs, lhs.y() * rhs, lhs.z() * rhs);
    }
    //=============================================================================
    Vec3 operator*(real_t lhs, const Vec3 &rhs)
    {
        return rhs * lhs;
    }
    //=============================================================================
    Vec3 operator/(const Vec3 &lhs, real_t rhs)
    {
        return lhs * (1 / rhs);
    }
    //=============================================================================
    Vec3 operator+(const Vec3 &lhs, real_t rhs)
    {
        return Vec3(lhs.x() + rhs, lhs.y() + rhs, lhs.z() + rhs);
    }
    //=============================================================================
    Vec3 operator+(real_t lhs, const Vec3 &rhs)
    {
        return rhs + lhs;
    }
    //=============================================================================
    Vec3 operator-(const Vec3 &lhs, real_t rhs)
    {
        return lhs + (-rhs);
    }
    //=============================================================================
    Vec3 operator-(real_t lhs, const Vec3 &rhs)
    {
        return rhs - lhs;
    }
    //=============================================================================
    Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs)
    {
        return Vec3(lhs.x() + rhs.x(), lhs.y() + rhs.y(), lhs.z() + rhs.z());
    }
    //=============================================================================
    Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs)
    {
        return Vec3(lhs.x() - rhs.x(), lhs.y() - rhs.y(), lhs.z() - rhs.z());
    }
}