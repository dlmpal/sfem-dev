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
    FieldCoefficient::FieldCoefficient(std::shared_ptr<Field> phi)
        : phi_(phi)
    {
    }
    //=============================================================================
    FieldCoefficient::FieldCoefficient(std::shared_ptr<const IndexMap> index_map,
                                       const std::vector<std::string> &components)
        : phi_(std::make_shared<Field>(index_map, components))
    {
    }
    //=============================================================================
    std::shared_ptr<Field> FieldCoefficient::field()
    {
        return phi_;
    }
    //=============================================================================
    std::shared_ptr<const Field> FieldCoefficient::field() const
    {
        return phi_;
    }
    //=============================================================================
    real_t &FieldCoefficient::operator()(int idx, int comp)
    {
        return (*phi_)(idx, comp);
    }
    //=============================================================================
    real_t FieldCoefficient::operator()(int idx, int comp) const
    {
        return (*phi_)(idx, comp);
    }
}