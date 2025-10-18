#include "coefficient.hpp"

namespace sfem
{
    //=============================================================================
    ConstantCoefficient::ConstantCoefficient(const std::vector<real_t> &value)
        : value_(value)
    {
    }
    //=============================================================================
    ConstantCoefficient::ConstantCoefficient(real_t value)
        : ConstantCoefficient(std::vector<real_t>{value})
    {
    }
    //=============================================================================
    real_t &ConstantCoefficient::operator()([[maybe_unused]] int, int comp)
    {
        return value_[comp];
    }
    //=============================================================================
    real_t ConstantCoefficient::operator()([[maybe_unused]] int, int comp) const
    {
        return value_[comp];
    }
    //=============================================================================
    FieldCoefficient::FieldCoefficient(Field &&field)
        : field_(std::move(field))
    {
    }
    //=============================================================================
    FieldCoefficient::FieldCoefficient(std::shared_ptr<IndexMap> index_map,
                                       const std::vector<std::string> &components)
        : field_(index_map, components)
    {
    }
    //=============================================================================
    Field &FieldCoefficient::field()
    {
        return field_;
    }
    //=============================================================================
    const Field &FieldCoefficient::field() const
    {
        return field_;
    }
    //=============================================================================
    real_t &FieldCoefficient::operator()(int idx, int comp)
    {
        return field_(idx, comp);
    }
    //=============================================================================
    real_t FieldCoefficient::operator()(int idx, int comp) const
    {
        return field_(idx, comp);
    }
}