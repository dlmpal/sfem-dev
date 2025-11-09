#pragma once

#include "fv_space.hpp"
#include "../field.hpp"

namespace sfem::fvm
{
    /// @brief Field defined on a finite volume space
    class FVField : public Field
    {
    public:
        /// @brief Create a FVField for a given space
        /// @param fv_space Finite volume space
        /// @param components Component names
        FVField(std::shared_ptr<const FVSpace> fv_space,
                const std::vector<std::string> &components);

        /// @brief Get the finite volume space
        std::shared_ptr<const FVSpace> space() const;

    private:
        /// @brief Finite volume space
        std::shared_ptr<const FVSpace> fv_space_;
    };

    /// @brief Create a finite volume field
    std::shared_ptr<FVField> create_field(std::shared_ptr<const FVSpace> V,
                                          const std::vector<std::string> &compoments);

    using FieldFunction = std::function<void(const std::array<real_t, 3> &pt,
                                             std::span<real_t> values, real_t time)>;

    /// @brief Explicitly evaluate a finite volume field
    /// @param phi Finite volume field
    /// @param func Evaluation function
    /// @param time Current simulation time
    void eval_field(FVField &phi, FieldFunction func, real_t time);
}