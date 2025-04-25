#include "xml.hpp"
#include "xml_utils.hpp"
#include <fstream>

namespace sfem::io::vtk::xml
{
    //=============================================================================
    static std::string real_t_str()
    {
        if constexpr (std::is_same_v<real_t, float>)
        {
            return "Float32";
        }
        else
        {
            return "Float64";
        }
    }
    //=============================================================================
    void write_vtu(const std::filesystem::path &filename,
                   const std::vector<int> &cell_types,
                   const graph::Connectivity &cell_to_node,
                   const std::vector<std::array<real_t, 3>> &points,
                   const std::vector<std::vector<real_t>> &cell_values,
                   const std::vector<std::string> &cell_names,
                   const std::vector<std::vector<real_t>> &node_values,
                   const std::vector<std::string> &node_names)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Indentation level
        int indent = 0;

        // Root element
        Element root(file, "VTKFile",
                     {{"type", "UnstructuredGrid"},
                      {"version", "0.1"},
                      {"byte_order", "BigEndian"}},
                     indent++);

        // File type element
        Element type(file, "UnstructuredGrid", {}, indent++);

        // Piece element
        Element piece(file, "Piece",
                      {{"NumberOfPoints", std::to_string(cell_to_node.n_secondary())},
                       {"NumberOfCells", std::to_string(cell_to_node.n_primary())}},
                      indent++);

        // Points
        {
            Element points_elem(file, "Points", {}, indent++);
            Element data_array(file, "DataArray",
                               {{"type", real_t_str()},
                                {"NumberOfComponents", "3"},
                                {"Format", "ascii"}},
                               indent++);
            for (const auto &point : points)
            {
                file << indent_string(std::format("{} {} {}\n", point[0], point[1], point[2]), indent);
            }
            indent -= 2;
        }

        // Cells
        {
            Element cells_elem(file, "Cells", {}, indent++);

            // Cell-to-node connectivity array
            {
                Element data_array(file, "DataArray",
                                   {{"type", "Int32"},
                                    {"Name", "connectivity"},
                                    {"Format", "ascii"}},
                                   indent++);
                for (int i = 0; i < cell_to_node.n_primary(); i++)
                {
                    file << indent_string("", indent);
                    for (int j : cell_to_node.links(i))
                    {
                        file << j << " ";
                    }
                    file << "\n";
                }
                indent--;
            }

            // Cell-to-node connectivity offsets
            {
                Element data_array(file, "DataArray",
                                   {{"type", "Int32"},
                                    {"Name", "offsets"},
                                    {"Format", "ascii"}},
                                   indent++);
                for (int i = 0; i < cell_to_node.n_primary(); i++)
                {
                    int offset = cell_to_node.offset(i) + cell_to_node.n_links(i);
                    file << indent_string(std::to_string(offset), indent) << "\n";
                }
                indent--;
            }

            // Cell types
            {
                Element data_array(file, "DataArray",
                                   {{"type", "Int32"},
                                    {"Name", "types"},
                                    {"Format", "ascii"}},
                                   indent++);
                for (int i = 0; i < cell_to_node.n_primary(); i++)
                {
                    file << indent_string(std::to_string(cell_types[i]), indent) << "\n";
                }
                indent--;
            }
            indent--;
        }

        // Cell data
        {
            Element cell_data(file, "CellData", {}, indent++);
            for (std::size_t i = 0; i < cell_values.size(); i++)
            {
                Element data_array(file, "DataArray",
                                   {{"type", real_t_str()},
                                    {"Name", cell_names[i]}},
                                   indent++);
                for (std::size_t j = 0; j < cell_values[i].size(); j++)
                {
                    file << indent_string(std::format("{}\n", cell_values[i][j]), indent);
                }
                indent--;
            }
            indent--;
        }

        // Point data
        {
            Element point_data(file, "PointData", {}, indent++);
            for (std::size_t i = 0; i < node_values.size(); i++)
            {
                Element data_array(file, "DataArray",
                                   {{"type", real_t_str()},
                                    {"Name", node_names[i]}},
                                   indent++);
                for (std::size_t j = 0; j < node_values[i].size(); j++)
                {
                    file << indent_string(std::format("{}\n", node_values[i][j]), indent);
                }
                indent--;
            }
            indent--;
        }
    }
    //=============================================================================
    void write_pvtu(const std::filesystem::path &filename,
                    const std::vector<std::filesystem::path> &sources,
                    const std::vector<std::string> &cell_names,
                    const std::vector<std::string> &node_names)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Indentation level
        int indent = 0;

        // Root element
        Element root(file, "VTKFile",
                     {{"type", "PUnstructuredGrid"},
                      {"version", "0.1"},
                      {"byte_order", "BigEndian"}},
                     indent++);

        // File type element
        Element type(file, "PUnstructuredGrid",
                     {{"GhostLevel", "0"}}, indent++);

        // Points
        {
            Element points(file, "PPoints", {}, indent++);
            Element data_array(file, "PDataArray",
                               {{"type", real_t_str()},
                                {"NumberOfComponents", "3"},
                                {"Format", "ascii"}},
                               indent--);
        }

        // Cells
        {
            Element cells(file, "PCells", {}, indent++);

            // Cell-to-node connectivity array
            {
                Element data_array(file, "PDataArray",
                                   {{"type", "Int32"},
                                    {"Name", "connectivity"},
                                    {"Format", "ascii"}},
                                   indent);
            }

            // Cell-to-node connectivity offsets
            {
                Element data_array(file, "PDataArray",
                                   {{"type", "Int32"},
                                    {"Name", "offsets"},
                                    {"Format", "ascii"}},
                                   indent);
            }

            // Cell types
            {
                Element data_array(file, "PDataArray",
                                   {{"type", "Int32"},
                                    {"Name", "types"},
                                    {"Format", "ascii"}},
                                   indent--);
            }
        }

        // Cell data
        {
            Element cell_data(file, "PCellData", {}, indent++);
            for (std::size_t i = 0; i < cell_names.size(); i++)
            {
                auto elem = create_empty_tag("PDataArray", {{"type", real_t_str()},
                                                            {"Name", cell_names[i]}});
                file << indent_string(elem, indent) << "\n";
            }
            indent--;
        }

        // Point data
        {
            Element cell_data(file, "PPointData", {}, indent++);
            for (std::size_t i = 0; i < node_names.size(); i++)
            {
                auto elem = create_empty_tag("PDataArray", {{"type", real_t_str()},
                                                            {"Name", node_names[i]}});
                file << indent_string(elem, indent) << "\n";
            }
            indent--;
        }

        // Sources
        {
            for (const auto &source : sources)
            {
                auto elem = create_empty_tag("Piece", {{"Source", source.string()}});
                file << indent_string(elem, indent) << "\n";
            }
        }
    }
}