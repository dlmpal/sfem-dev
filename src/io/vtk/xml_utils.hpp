#pragma once

#include <array>
#include <vector>
#include <string>
#include <fstream>

namespace sfem::io::vtk::xml
{
    /// @brief Start tag for an element
    std::string create_start_tag(const std::string &tag,
                                 const std::vector<std::array<std::string, 2>> &attributes = {});

    /// @brief End tag for an element
    std::string create_end_tag(const std::string &tag);

    /// @brief Create a tag for an empty element
    std::string create_empty_tag(const std::string &tag,
                                 const std::vector<std::array<std::string, 2>> &attributes = {});

    /// @brief Indent a string by the given number of tabs
    std::string indent_string(const std::string &str, int indent);

    class Element
    {
    public:
        Element(std::ofstream &file, const std::string &tag,
                const std::vector<std::array<std::string, 2>> &attributes = {},
                int indent = 0);
        ~Element();

    private:
        std::ofstream *file_;
        std::string tag_;
        int indent_;
    };
}