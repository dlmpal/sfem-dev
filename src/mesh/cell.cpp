#include "cell.hpp"
#include <map>

namespace sfem::mesh
{
    //=============================================================================
    std::string cell_type_str(CellType cell_type)
    {
        switch (cell_type)
        {
        case CellType::point:
            return "point";
        case CellType::line:
            return "line";
        case CellType::triangle:
            return "triangle";
        case CellType::quadrilateral:
            return "quadrilateral";
        case CellType::tetrahedron:
            return "tetrahedron";
        case CellType::hexahedron:
            return "hexahedron";
        case CellType::prism:
            return "prism";
        default:
            return "invalid";
        }
    }
    //=============================================================================
    int cell_dim(CellType cell_type)
    {
        switch (cell_type)
        {
        case CellType::point:
            return 0;
        case CellType::line:
            return 1;
        case CellType::triangle:
            return 2;
        case CellType::quadrilateral:
            return 2;
        case CellType::tetrahedron:
            return 3;
        case CellType::hexahedron:
            return 3;
        case CellType::prism:
            return 3;
        default:
            SFEM_BAD_CELL_ERROR(cell_type);
            return -1;
        }
    }
    //=============================================================================
    int cell_num_nodes(CellType cell_type)
    {
        switch (cell_type)
        {
        case CellType::point:
            return 1;
        case CellType::line:
            return 2;
        case CellType::triangle:
            return 3;
        case CellType::quadrilateral:
            return 4;
        case CellType::tetrahedron:
            return 4;
        case CellType::hexahedron:
            return 8;
        case CellType::prism:
            return 6;
        default:
            SFEM_BAD_CELL_ERROR(cell_type);
            return -1;
        }
    }
    //=============================================================================
    int cell_num_edges(CellType cell_type)
    {
        switch (cell_type)
        {
        case CellType::point:
            return 0;
        case CellType::line:
            return 1;
        case CellType::triangle:
            return 3;
        case CellType::quadrilateral:
            return 4;
        case CellType::tetrahedron:
            return 6;
        case CellType::hexahedron:
            return 12;
        case CellType::prism:
            return 9;
        default:
            SFEM_BAD_CELL_ERROR(cell_type);
            return -1;
        }
    }
    //=============================================================================
    int cell_num_faces(CellType cell_type)
    {
        switch (cell_type)
        {
        case CellType::point:
            return 0;
        case CellType::line:
            return 0;
        case CellType::triangle:
            return 1;
        case CellType::quadrilateral:
            return 1;
        case CellType::tetrahedron:
            return 4;
        case CellType::hexahedron:
            return 6;
        case CellType::prism:
            return 5;
        default:
            SFEM_BAD_CELL_ERROR(cell_type);
            return -1;
        }
    }
    //=============================================================================
    std::array<int, 2> cell_edge_ordering(CellType cell_type, int edge_idx)
    {
        SFEM_CHECK_INDEX(edge_idx, cell_num_edges(cell_type));

        static const std::map<std::pair<CellType, int>, std::array<int, 2>> orderings =
            {
                {{CellType::line, 0}, {0, 1}},

                {{CellType::triangle, 0}, {0, 1}},
                {{CellType::triangle, 1}, {1, 2}},
                {{CellType::triangle, 2}, {2, 0}},

                {{CellType::quadrilateral, 0}, {0, 1}},
                {{CellType::quadrilateral, 1}, {1, 2}},
                {{CellType::quadrilateral, 2}, {2, 3}},
                {{CellType::quadrilateral, 3}, {3, 0}},

                {{CellType::tetrahedron, 0}, {0, 1}},
                {{CellType::tetrahedron, 1}, {1, 2}},
                {{CellType::tetrahedron, 2}, {2, 0}},
                {{CellType::tetrahedron, 3}, {0, 3}},
                {{CellType::tetrahedron, 4}, {3, 2}},
                {{CellType::tetrahedron, 5}, {3, 1}},

                {{CellType::hexahedron, 0}, {0, 1}},
                {{CellType::hexahedron, 1}, {0, 3}},
                {{CellType::hexahedron, 2}, {0, 4}},
                {{CellType::hexahedron, 3}, {1, 2}},
                {{CellType::hexahedron, 4}, {1, 5}},
                {{CellType::hexahedron, 5}, {2, 3}},
                {{CellType::hexahedron, 6}, {2, 6}},
                {{CellType::hexahedron, 7}, {3, 7}},
                {{CellType::hexahedron, 8}, {4, 5}},
                {{CellType::hexahedron, 9}, {4, 7}},
                {{CellType::hexahedron, 10}, {5, 6}},
                {{CellType::hexahedron, 11}, {6, 7}},
            };

        return orderings.at({cell_type, edge_idx});
    }
    //=============================================================================
    CellType cell_face_type(CellType cell_type, int face_idx)
    {
        SFEM_CHECK_INDEX(face_idx, cell_num_faces(cell_type));

        static const std::map<std::pair<CellType, int>, CellType> face_type =
            {
                {{CellType::triangle, 0}, CellType::triangle},

                {{CellType::quadrilateral, 0}, CellType::quadrilateral},

                {{CellType::tetrahedron, 0}, CellType::triangle},
                {{CellType::tetrahedron, 1}, CellType::triangle},
                {{CellType::tetrahedron, 2}, CellType::triangle},
                {{CellType::tetrahedron, 3}, CellType::triangle},

                {{CellType::hexahedron, 0}, CellType::quadrilateral},
                {{CellType::hexahedron, 1}, CellType::quadrilateral},
                {{CellType::hexahedron, 2}, CellType::quadrilateral},
                {{CellType::hexahedron, 3}, CellType::quadrilateral},
                {{CellType::hexahedron, 4}, CellType::quadrilateral},
                {{CellType::hexahedron, 5}, CellType::quadrilateral},
            };

        return face_type.at({cell_type, face_idx});
    }
    //=============================================================================
    std::array<int, 4> cell_face_ordering(CellType cell_type, int face_idx)
    {
        SFEM_CHECK_INDEX(face_idx, cell_num_faces(cell_type));

        static const std::map<std::pair<CellType, int>, std::array<int, 4>> orderings =
            {
                {{CellType::triangle, 0}, {0, 1, 2}},
                {{CellType::quadrilateral, 0}, {0, 1, 2, 3}},

                {{CellType::tetrahedron, 0}, {1, 2, 3}},
                {{CellType::tetrahedron, 1}, {0, 2, 3}},
                {{CellType::tetrahedron, 2}, {0, 1, 3}},
                {{CellType::tetrahedron, 3}, {0, 1, 2}},

                {{CellType::hexahedron, 0}, {0, 3, 2, 1}},
                {{CellType::hexahedron, 1}, {0, 1, 5, 4}},
                {{CellType::hexahedron, 2}, {0, 4, 7, 3}},
                {{CellType::hexahedron, 3}, {1, 2, 6, 5}},
                {{CellType::hexahedron, 4}, {3, 7, 6, 2}},
                {{CellType::hexahedron, 5}, {4, 5, 6, 7}},
            };

        return orderings.at({cell_type, face_idx});
    }
}