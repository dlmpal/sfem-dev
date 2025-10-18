#include "xml.hpp"
#include "xml_utils.hpp"
#include "utils.hpp"
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
                   const std::vector<std::shared_ptr<const Field>> &cell_fields,
                   const std::vector<std::shared_ptr<const Field>> &node_fields)
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
                    std::vector<int> cell_nodes;
                    cell_nodes.assign(cell_to_node.links(i).cbegin(),
                                      cell_to_node.links(i).cend());
                    cell_node_ordering_to_vtk(cell_types[i], cell_nodes);
                    for (int j = 0; j < cell_to_node.n_links(j); j++)
                    {
                        file << cell_nodes[j] << " ";
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
            for (const auto &field : cell_fields)
            {
                for (const auto &comp : field->components())
                {
                    Element data_array_elem(file, "DataArray",
                                            {{"type", real_t_str()},
                                             {"Name", comp}},
                                            indent++);
                    const int comp_idx = field->comp_idx(comp);
                    for (int i = 0; i < field->n_local(); i++)
                    {
                        file << indent_string(std::format("{}\n", (*field)(i, comp_idx)), indent);
                    }
                    indent--;
                }
            }
            indent--;
        }

        // Point data
        {
            Element point_data_elem(file, "PointData", {}, indent++);
            for (const auto &field : node_fields)
            {
                for (const auto &comp : field->components())
                {
                    Element data_array_elem(file, "DataArray",
                                            {{"type", real_t_str()},
                                             {"Name", comp}},
                                            indent++);
                    const int comp_idx = field->comp_idx(comp);
                    for (int i = 0; i < field->n_local(); i++)
                    {
                        file << indent_string(std::format("{}\n", (*field)(i, comp_idx)), indent);
                    }
                    indent--;
                }
            }
            indent--;
        }
    }
    //=============================================================================
    void write_pvtu(const std::filesystem::path &filename,
                    const std::vector<std::filesystem::path> &sources,
                    const std::vector<std::shared_ptr<const Field>> &cell_fields,
                    const std::vector<std::shared_ptr<const Field>> &node_fields)
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
            for (const auto &field : cell_fields)
            {
                for (const auto &comp : field->components())
                {
                    auto elem = create_empty_tag("PDataArray", {{"type", real_t_str()},
                                                                {"Name", comp}});
                    file << indent_string(elem, indent) << "\n";
                }
            }
            indent--;
        }

        // Point data
        {
            Element point_data_elem(file, "PPointData", {}, indent++);
            for (const auto &field : node_fields)
            {
                for (const auto &comp : field->components())
                {
                    auto elem = create_empty_tag("PDataArray", {{"type", real_t_str()},
                                                                {"Name", comp}});
                    file << indent_string(elem, indent) << "\n";
                }
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