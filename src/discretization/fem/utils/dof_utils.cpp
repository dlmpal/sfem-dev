#include "dof_utils.hpp"
#include "../../../geo/utils.hpp"
#include "../../../base/error.hpp"

namespace sfem::fem::dof
{
    //=============================================================================
    int cell_num_dof(mesh::CellType cell_type, int order)
    {
        if (order == 0)
        {
            return 1;
        }

        switch (cell_type)
        {
        case mesh::CellType::point:
            return 1;
        case mesh::CellType::line:
            return order + 1;
        case mesh::CellType::triangle:
            return (order + 1) * (order + 2) / 2;
        case mesh::CellType::quadrilateral:
            return (order + 1) * (order + 1);
        case mesh::CellType::tetrahedron:
            return (order + 1) * (order + 2) * (order + 3) / 6;
        case mesh::CellType::hexahedron:
            return (order + 1) * (order + 1) * (order + 1);
        case mesh::CellType::prism:
            return (order + 1) * (order + 1) * (order + 2) / 2;
        default:
            SFEM_BAD_CELL_ERROR(cell_type);
            return 0;
        }
    }
    //=============================================================================
    int cell_num_internal_dof(mesh::CellType cell_type, int order)
    {
        if (order == 0)
        {
            return 1;
        }

        switch (cell_type)
        {
        case mesh::CellType::point:
            return 1;
        case mesh::CellType::line:
            return order - 1;
        case mesh::CellType::triangle:
            return (order - 1) * (order - 2) / 2;
        case mesh::CellType::quadrilateral:
            return (order - 1) * (order - 1);
        case mesh::CellType::tetrahedron:
            return (order - 1) * (order - 2) * (order - 3) / 6;
        case mesh::CellType::hexahedron:
            return (order - 1) * (order - 1) * (order - 1);
        case mesh::CellType::prism:
            return (order - 1) * (order - 1) * (order - 2) / 2;
        default:
            SFEM_BAD_CELL_ERROR(cell_type);
            return 0;
        }
    }
    //=============================================================================
    void compute_line_dof_points(int order, std::vector<std::array<real_t, 3>> &points)
    {
        // Intneral DoF
        for (int i = 0; i < cell_num_internal_dof(mesh::CellType::line, order); i++)
        {
            points.emplace_back(geo::compute_line_nth_point(static_cast<int>(i) + 1,
                                                            order,
                                                            points[0],
                                                            points[1]));
        }
    }
    //=============================================================================
    void compute_triangle_dof_points(int order, std::vector<std::array<real_t, 3>> &points)
    {
        // No. corner nodes, edges and dof per edge
        int n_nodes = cell_num_nodes(mesh::CellType::triangle);
        int n_edges = cell_num_edges(mesh::CellType::triangle);
        int n_dof_edge = dof::cell_num_internal_dof(mesh::CellType::line, order);

        // Edge DoF
        for (int i = 0; i < n_edges; i++)
        {
            auto edge_ordering = mesh::cell_edge_ordering(mesh::CellType::triangle, i);
            auto edge_point_1 = points[edge_ordering[0]];
            auto edge_point_2 = points[edge_ordering[1]];

            for (int j = 0; j < n_dof_edge; j++)
            {
                points.emplace_back(geo::compute_line_nth_point(j + 1,
                                                                order,
                                                                edge_point_1,
                                                                edge_point_2));
            }
        }

        // Quick access to edge points
        std::span<std::array<real_t, 3>> edge_points(points.begin() + n_nodes,
                                                     n_edges * n_dof_edge);

        // Internal DoF
        for (int i = 0; i < n_dof_edge - 1; i++)
        {
            auto edge_1_point_1 = edge_points[i];
            auto edge_1_point_2 = edge_points[n_dof_edge * 2 - i - 1];

            for (int j = 0; j < n_dof_edge - i - 1; j++)
            {
                auto edge_2_point_1 = edge_points[n_dof_edge * 3 - j - 1];
                auto edge_2_point_2 = edge_points[n_dof_edge + j];

                points.emplace_back(geo::compute_line_intersection(edge_1_point_1,
                                                                   edge_1_point_2,
                                                                   edge_2_point_1,
                                                                   edge_2_point_2));
            }
        }
    }
    //=============================================================================
    void compute_quad_dof_points(int order, std::vector<std::array<real_t, 3>> &points)
    {
        // No. corner nodes, edges and dof per edge
        int n_nodes = mesh::cell_num_nodes(mesh::CellType::quadrilateral);
        int n_edges = mesh::cell_num_edges(mesh::CellType::quadrilateral);
        int n_dof_edge = dof::cell_num_internal_dof(mesh::CellType::line, order);

        // Edge DoF
        for (int i = 0; i < n_edges; i++)
        {
            auto edge_ordering = mesh::cell_edge_ordering(mesh::CellType::quadrilateral, i);
            auto edge_point_1 = points[edge_ordering[0]];
            auto edge_point_2 = points[edge_ordering[1]];

            for (int j = 0; j < n_dof_edge; j++)
            {
                points.emplace_back(geo::compute_line_nth_point(j + 1,
                                                                order,
                                                                edge_point_1,
                                                                edge_point_2));
            }
        }

        // Quick access to edge points
        std::span<std::array<real_t, 3>> edge_points(points.begin() + n_nodes,
                                                     n_edges * n_dof_edge);

        // Internal DoF
        for (int i = 0; i < n_dof_edge; i++)
        {
            auto edge_1_point_1 = edge_points[n_dof_edge + i];
            auto edge_1_point_2 = edge_points[n_dof_edge * n_edges - i - 1];

            for (int j = 0; j < n_dof_edge; j++)
            {
                auto edge_2_point_1 = edge_points[j];
                auto edge_2_point_2 = edge_points[n_dof_edge * (n_edges - 1) - j - 1];

                points.emplace_back(geo::compute_line_intersection(edge_1_point_1,
                                                                   edge_1_point_2,
                                                                   edge_2_point_1,
                                                                   edge_2_point_2));
            }
        }
    }
    //=============================================================================
    void compute_tet_dof_points(int order, std::vector<std::array<real_t, 3>> &points)
    {
        // No. corner nodes, edges and dof per edge
        int n_nodes = cell_num_nodes(mesh::CellType::tetrahedron);
        int n_edges = cell_num_edges(mesh::CellType::tetrahedron);
        int n_dof_edge = dof::cell_num_internal_dof(mesh::CellType::line, order);

        // Edge DoF
        for (int i = 0; i < n_edges; i++)
        {
            auto edge_ordering = mesh::cell_edge_ordering(mesh::CellType::tetrahedron, i);
            auto edge_point_1 = points[edge_ordering[0]];
            auto edge_point_2 = points[edge_ordering[1]];

            for (int j = 0; j < n_dof_edge; j++)
            {
                points.emplace_back(geo::compute_line_nth_point(j + 1,
                                                                order,
                                                                edge_point_1,
                                                                edge_point_2));
            }
        }

        // Quick access to edge points
        std::span<std::array<real_t, 3>> edge_points(points.begin() + n_nodes,
                                                     n_edges * n_dof_edge);

        /// @todo Face DoF
        /// @todo Internal DoF
    }
    //=============================================================================
    void compute_hex_dof_points(int order, std::vector<std::array<real_t, 3>> &points)
    {
        // No. corner nodes, edges and dof per edge
        int n_nodes = cell_num_nodes(mesh::CellType::hexahedron);
        int n_edges = cell_num_edges(mesh::CellType::hexahedron);
        int n_dof_edge = dof::cell_num_internal_dof(mesh::CellType::line, order);

        // Edge DoF
        for (int i = 0; i < n_edges; i++)
        {
            auto edge_ordering = mesh::cell_edge_ordering(mesh::CellType::hexahedron, i);
            auto edge_point_1 = points[edge_ordering[0]];
            auto edge_point_2 = points[edge_ordering[1]];

            for (int j = 0; j < n_dof_edge; j++)
            {
                points.emplace_back(geo::compute_line_nth_point(j + 1,
                                                                order,
                                                                edge_point_1,
                                                                edge_point_2));
            }
        }

        // Quick access to edge points
        std::span<std::array<real_t, 3>> edge_points(points.begin() + n_nodes,
                                                     n_edges * n_dof_edge);

        /// @todo Face DoF
        /// @todo Internal DoF
    }
    //=============================================================================
    void compute_cell_dof_points(mesh::CellType cell_type, int order,
                                 std::vector<std::array<real_t, 3>> &points)
    {
        // Return early for constant/linear elements
        if (order <= 1)
        {
            return;
        }

        switch (cell_type)
        {
        case mesh::CellType::line:
            compute_line_dof_points(order, points);
            break;
        case mesh::CellType::triangle:
            compute_triangle_dof_points(order, points);
            break;
        case mesh::CellType::quadrilateral:
            compute_quad_dof_points(order, points);
            break;
        case mesh::CellType::tetrahedron:
            compute_tet_dof_points(order, points);
            break;
        case mesh::CellType::hexahedron:
            compute_tet_dof_points(order, points);
            break;
        default:
            SFEM_BAD_CELL_ERROR(cell_type);
        }
    }
}