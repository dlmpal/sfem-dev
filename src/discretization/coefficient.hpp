#pragma once

#include "field.hpp"

namespace sfem
{
    /// @brief Coefficients represent various physical quantities
    /// present in the underlying PDE. They are used by both FEM and FVM
    /// assembly kernels
    class Coefficient
    {
    public:
        virtual real_t &operator()(int idx, int comp) = 0;
        virtual real_t operator()(int idx, int comp) const = 0;
        virtual int n_comp() const = 0;
    };

    /// @brief A coefficient thats constant everywhere
    class ConstantCoefficient : public Coefficient
    {
    public:
        ConstantCoefficient(const std::vector<real_t> &value);
        ConstantCoefficient(real_t value);
        real_t &operator()(int idx, int comp) override;
        real_t operator()(int idx, int comp) const override;
        int n_comp() const override;

    private:
        /// @brief Constant value
        std::vector<real_t> value_;
    };

    /// @brief A coefficient whose values are defined by an underlying field
    class FieldCoefficient : public Coefficient
    {
    public:
        FieldCoefficient(std::shared_ptr<Field> phi);
        FieldCoefficient(std::shared_ptr<const IndexMap>,
                         const std::vector<std::string> &components);

        std::shared_ptr<Field> field();
        std::shared_ptr<const Field> field() const;

        real_t &operator()(int idx, int comp) override;
        real_t operator()(int idx, int comp) const override;
        int n_comp() const override;

    protected:
        std::shared_ptr<Field> phi_;
    };
}