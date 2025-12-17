#pragma once

#include "../base/config.hpp"
#include <functional>
#include <string>

namespace sfem::io
{
    class IOField
    {
    public:
        using ValueAccessor = std::function<real_t(int idx, int comp)>;

        IOField(const std::vector<std::string> &components, ValueAccessor func);

        std::vector<std::string> components() const;
        int comp_idx(const std::string &component) const;
        real_t operator()(int idx, int comp) const;

    private:
        std::vector<std::string> components_;
        ValueAccessor values_;
    };
}