#pragma once

#include "../elements/fe.hpp"
#include "../../coefficient.hpp"
#include "../fe_function.hpp"

namespace sfem::fem::kernels
{
    class LinearElasticity2D
    {
    public:
        LinearElasticity2D(std::shared_ptr<const Coefficient> E,
                           std::shared_ptr<const Coefficient> nu,
                           std::shared_ptr<const Coefficient> thick);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data) const;

    private:
        std::shared_ptr<const Coefficient> E_;
        std::shared_ptr<const Coefficient> nu_;
        std::shared_ptr<const Coefficient> thick_;
    };

    class InertialLoad2D
    {
    public:
        InertialLoad2D(std::shared_ptr<const Coefficient> thick,
                       std::shared_ptr<const Coefficient> rho,
                       const std::array<real_t, 2> &g);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data);

    private:
        std::shared_ptr<const Coefficient> thick_;
        std::shared_ptr<const Coefficient> rho_;
        std::array<real_t, 2> g_;
    };

    class PressureLoad2D
    {
    public:
        PressureLoad2D(std::shared_ptr<const Coefficient> thick,
                       std::shared_ptr<const Coefficient> pressure);

        la::DenseMatrix operator()(int facet_idx,
                                   const fem::FEData &data,
                                   const geo::Vec3 &normal);

    private:
        std::shared_ptr<const Coefficient> thick_;
        std::shared_ptr<const Coefficient> pressure_;
    };

    class LinearElasticity3D
    {
    public:
        LinearElasticity3D(std::shared_ptr<const Coefficient> E,
                           std::shared_ptr<const Coefficient> nu);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data) const;

    private:
        std::shared_ptr<const Coefficient> E_;
        std::shared_ptr<const Coefficient> nu_;
    };

    class VonMises3D
    {
    public:
        VonMises3D(std::shared_ptr<const Coefficient> E,
                   std::shared_ptr<const Coefficient> nu,
                   std::shared_ptr<const FEField> U);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data) const;

    private:
        std::shared_ptr<const Coefficient> E_;
        std::shared_ptr<const Coefficient> nu_;
        std::shared_ptr<const FEField> U_;
    };

    class InertialLoad3D
    {
    public:
        InertialLoad3D(std::shared_ptr<const Coefficient> rho, const std::array<real_t, 3> &g);
        la::DenseMatrix operator()(int cell_idx, const fem::FEData &data);

    private:
        std::shared_ptr<const Coefficient> rho_;
        std::array<real_t, 3> g_;
    };

    class PressureLoad3D
    {
    public:
        PressureLoad3D(std::shared_ptr<const Coefficient> pressure);

        la::DenseMatrix operator()(int facet_idx,
                                   const fem::FEData &data,
                                   const geo::Vec3 &normal);

    private:
        std::shared_ptr<const Coefficient> pressure_;
    };
}