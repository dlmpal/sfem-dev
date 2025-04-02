#pragma once

#include <string>

namespace sfem::mesh
{
    class Region
    {
    public:
        /// @brief Create a region
        Region(const std::string &name, int tag, int dim);

        /// @brief Get the region's name
        std::string name() const;

        /// @brief Get the region's tag
        int tag() const;

        /// @brief Get the region's dimension
        int dim() const;

    private:
        /// @brief Name
        std::string name_;

        /// @brief Integer tag
        int tag_;

        /// @brief Physical dimension
        int dim_;
    };
}
