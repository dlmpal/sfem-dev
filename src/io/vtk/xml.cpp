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
                   const std::vector<std::pair<std::string, std::span<real_t>>> &cell_data,
                   const std::vector<std::pair<std::string, std::span<real_t>>> &node_data)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Indentation level
        int indent = 0;

        // Root element
        Element root_elem(file, "VTKFile",
                          {{"type", "UnstructuredGrid"},
                           {"version", "0.1"},
                           {"byte_order", "BigEndian"}},
                          indent++);

        // File type element
        Element type_elem(file, "UnstructuredGrid", {}, indent++);

        // Piece element
        Element piece_elem(file, "Piece",
                           {{"NumberOfPoints", std::to_string(cell_to_node.n_secondary())},
                            {"NumberOfCells", std::to_string(cell_to_node.n_primary())}},
                           indent++);

        // Points
        {
            Element points_elem(file, "Points", {}, indent++);
            Element data_array_elem(file, "DataArray",
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
                Element data_array_elem(file, "DataArray",
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
                Element data_array_elem(file, "DataArray",
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
                Element data_array_elem(file, "DataArray",
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
            Element cell_data_elem(file, "CellData", {}, indent++);
            for (const auto &[name, data] : cell_data)
            {
                Element data_array_elem(file, "DataArray",
                                        {{"type", real_t_str()},
                                         {"Name", name}},
                                        indent++);
                for (std::size_t i = 0; i < data.size(); i++)
                {
                    file << indent_string(std::format("{}\n", data[i]), indent);
                }
                indent--;
            }
            indent--;
        }

        // Point data
        {
            Element point_data_elem(file, "PointData", {}, indent++);
            for (const auto &[name, data] : node_data)
            {
                Element data_array_elem(file, "DataArray",
                                        {{"type", real_t_str()},
                                         {"Name", name}},
                                        indent++);
                for (std::size_t i = 0; i < data.size(); i++)
                {
                    file << indent_string(std::format("{}\n", data[i]), indent);
                }
                indent--;
            }
            indent--;
        }
    }
    //=============================================================================
    void write_pvtu(const std::filesystem::path &filename,
                    const std::vector<std::filesystem::path> &sources,
                    const std::vector<std::pair<std::string, std::span<real_t>>> &cell_data,
                    const std::vector<std::pair<std::string, std::span<real_t>>> &node_data)
    {
        std::ofstream file(filename);
        SFEM_CHECK_FILE_OPEN(file, filename);

        // Indentation level
        int indent = 0;

        // Root element
        Element root_elem(file, "VTKFile",
                          {{"type", "PUnstructuredGrid"},
                           {"version", "0.1"},
                           {"byte_order", "BigEndian"}},
                          indent++);

        // File type element
        Element type_elem(file, "PUnstructuredGrid",
                          {{"GhostLevel", "0"}}, indent++);

        // Points
        {
            Element points_elem(file, "PPoints", {}, indent++);
            Element data_array_elem(file, "PDataArray",
                                    {{"type", real_t_str()},
                                     {"NumberOfComponents", "3"},
                                     {"Format", "ascii"}},
                                    indent--);
        }

        // Cells
        {
            Element cells_elem(file, "PCells", {}, indent++);

            // Cell-to-node connectivity array
            {
                Element data_array_elem(file, "PDataArray",
                                        {{"type", "Int32"},
                                         {"Name", "connectivity"},
                                         {"Format", "ascii"}},
                                        indent);
            }

            // Cell-to-node connectivity offsets
            {
                Element data_array_elem(file, "PDataArray",
                                        {{"type", "Int32"},
                                         {"Name", "offsets"},
                                         {"Format", "ascii"}},
                                        indent);
            }

            // Cell types
            {
                Element data_array_elem(file, "PDataArray",
                                        {{"type", "Int32"},
                                         {"Name", "types"},
                                         {"Format", "ascii"}},
                                        indent--);
            }
        }

        // Cell data
        {
            Element cell_data_elem(file, "PCellData", {}, indent++);
            for (const auto &[name, _] : cell_data)
            {
                auto elem = create_empty_tag("PDataArray", {{"type", real_t_str()},
                                                            {"Name", name}});
                file << indent_string(elem, indent) << "\n";
            }
            indent--;
        }

        // Point data
        {
            Element point_data_elem(file, "PPointData", {}, indent++);
            for (const auto &[name, _] : node_data)
            {
                auto elem = create_empty_tag("PDataArray", {{"type", real_t_str()},
                                                            {"Name", name}});
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