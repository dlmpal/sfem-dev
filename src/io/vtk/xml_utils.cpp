#include "xml_utils.hpp"
#include <format>

namespace sfem::io::vtk::xml
{
    //=============================================================================
    std::string create_start_tag(const std::string &tag,
                                 const std::vector<std::array<std::string, 2>> &attributes)
    {
        // Concatenate all attributes in a single string
        std::string metadata = "";
        for (const auto &[key, value] : attributes)
        {
            metadata += std::format(" {}=\"{}\"", key, value);
        }

        return std::format("<{}{}>", tag, metadata);
    }
    //=============================================================================
    std::string create_end_tag(const std::string &tag)
    {
        return std::format("</{}>", tag);
    }
    //=============================================================================
    std::string create_empty_tag(const std::string &tag,
                                 const std::vector<std::array<std::string, 2>> &attributes)
    {
        // Concatenate all attributes in a single string
        std::string metadata = "";
        for (const auto &[key, value] : attributes)
        {
            metadata += std::format(" {}=\"{}\"", key, value);
        }

        return std::format("<{}{}/>", tag, metadata);
    }
    //=============================================================================
    std::string indent_string(const std::string &str, int indent_level)
    {
        return std::string(indent_level, '\t') + str;
    }
    //=============================================================================
    Element::Element(std::ofstream &file, const std::string &tag,
                     const std::vector<std::array<std::string, 2>> &attributes,
                     int indent)
        : file_(&file),
          tag_(tag),
          indent_(indent)
    {
        *file_ << indent_string(create_start_tag(tag_, attributes), indent_) << "\n";
    }
    //=============================================================================
    Element::~Element()
    {
        *file_ << indent_string(create_end_tag(tag_), indent_) << "\n";
    }
}