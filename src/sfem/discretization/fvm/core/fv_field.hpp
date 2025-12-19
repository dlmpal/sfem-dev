#pragma once

#include <sfem/discretization/fvm/core/fv_space.hpp>
#include <sfem/discretization/fvm/core/fv_bc.hpp>
#include <sfem/discretization/fvm/core/fv_gradient.hpp>

// Forward declaration
namespace sfem::la
{
    class Vector;
}

namespace sfem::fvm
{
    class IField
    {
    public:
        IField(const std::vector<std::string> &components);
        virtual ~IField() = default;

        std::vector<std::string> components() const;
        int n_comp() const;

        virtual real_t &cell_value(int cell_idx, int comp_idx = 0) = 0;
        virtual real_t cell_value(int cell_idx, int comp_idx = 0) const = 0;

        virtual real_t facet_value(int facet_idx, int comp_idx = 0) const = 0;

        virtual geo::Vec3 cell_grad(int cell_idx, int comp_idx = 0) const = 0;
        virtual geo::Vec3 facet_grad(int cell_idx, int comp_idx = 0) const = 0;

    protected:
        std::vector<std::string> components_;
    };

    class ConstantField : public IField
    {
    public:
        ConstantField(const std::vector<std::string> &components,
                      const std::vector<real_t> &value);

        ConstantField(const std::string &name, real_t value);

        real_t &cell_value(int cell_idx, int comp_idx = 0) override;
        real_t cell_value(int cell_idx, int comp_idx = 0) const override;

        real_t facet_value(int facet_idx, int comp_idx = 0) const override;

        geo::Vec3 cell_grad(int cell_idx, int comp_idx = 0) const override;
        geo::Vec3 facet_grad(int cell_idx, int comp_idx = 0) const override;

    private:
        std::vector<real_t> value_;
    };

    class FVField : public IField
    {
    public:
        FVField(std::shared_ptr<const FVSpace> V,
                const std::vector<std::string> &components,
                GradientMethod gradient_method = GradientMethod::none);

        std::shared_ptr<const FVSpace> space() const;

        FVBC &boundary_condition();
        const FVBC &boundary_condition() const;

        la::Vector &values();
        const la::Vector &values() const;

        GradientMethod grad_method() const;

        la::Vector &grad();
        const la::Vector &grad() const;

        real_t &cell_value(int cell_idx, int comp_idx = 0) override;
        real_t cell_value(int cell_idx, int comp_idx = 0) const override;

        real_t facet_value(int facet_idx, int comp_idx = 0) const override;

        geo::Vec3 cell_grad(int cell_idx, int comp_idx = 0) const override;
        geo::Vec3 facet_grad(int facet_idx, int comp_idx = 0) const override;

        void update_gradient();

    private:
        std::shared_ptr<const FVSpace> V_;

        std::shared_ptr<FVBC> bc_;

        std::shared_ptr<la::Vector> values_;

        GradientMethod gradient_method_;

        std::shared_ptr<la::Vector> grad_;
    };
}