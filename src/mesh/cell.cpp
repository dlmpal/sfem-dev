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
        case CellType::quad:
            return "quadrilateral";
        case CellType::tet:
            return "tetrahedron";
        case CellType::hex:
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
        case CellType::quad:
            return 2;
        case CellType::tet:
            return 3;
        case CellType::hex:
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
        case CellType::quad:
            return 4;
        case CellType::tet:
            return 4;
        case CellType::hex:
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
        case CellType::quad:
            return 4;
        case CellType::tet:
            return 6;
        case CellType::hex:
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
        case CellType::quad:
            return 1;
        case CellType::tet:
            return 4;
        case CellType::hex:
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

                {{CellType::quad, 0}, {0, 1}},
                {{CellType::quad, 1}, {1, 2}},
                {{CellType::quad, 2}, {2, 3}},
                {{CellType::quad, 3}, {3, 0}},

                {{CellType::tet, 0}, {0, 1}},
                {{CellType::tet, 1}, {1, 2}},
                {{CellType::tet, 2}, {2, 0}},
                {{CellType::tet, 3}, {0, 3}},
                {{CellType::tet, 4}, {3, 2}},
                {{CellType::tet, 5}, {3, 1}},

                {{CellType::hex, 0}, {0, 1}},
                {{CellType::hex, 1}, {0, 3}},
                {{CellType::hex, 2}, {0, 4}},
                {{CellType::hex, 3}, {1, 2}},
                {{CellType::hex, 4}, {1, 5}},
                {{CellType::hex, 5}, {2, 3}},
                {{CellType::hex, 6}, {2, 6}},
                {{CellType::hex, 7}, {3, 7}},
                {{CellType::hex, 8}, {4, 5}},
                {{CellType::hex, 9}, {4, 7}},
                {{CellType::hex, 10}, {5, 6}},
                {{CellType::hex, 11}, {6, 7}},
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

                {{CellType::quad, 0}, CellType::quad},

                {{CellType::tet, 0}, CellType::triangle},
                {{CellType::tet, 1}, CellType::triangle},
                {{CellType::tet, 2}, CellType::triangle},
                {{CellType::tet, 3}, CellType::triangle},

                {{CellType::hex, 0}, CellType::quad},
                {{CellType::hex, 1}, CellType::quad},
                {{CellType::hex, 2}, CellType::quad},
                {{CellType::hex, 3}, CellType::quad},
                {{CellType::hex, 4}, CellType::quad},
                {{CellType::hex, 5}, CellType::quad},
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
                {{CellType::quad, 0}, {0, 1, 2, 3}},

                {{CellType::tet, 0}, {1, 2, 3}},
                {{CellType::tet, 1}, {0, 2, 3}},
                {{CellType::tet, 2}, {0, 1, 3}},
                {{CellType::tet, 3}, {0, 1, 2}},

                {{CellType::hex, 0}, {0, 3, 2, 1}},
                {{CellType::hex, 1}, {0, 1, 5, 4}},
                {{CellType::hex, 2}, {0, 4, 7, 3}},
                {{CellType::hex, 3}, {1, 2, 6, 5}},
                {{CellType::hex, 4}, {3, 7, 6, 2}},
                {{CellType::hex, 5}, {4, 5, 6, 7}},
            };

        return orderings.at({cell_type, face_idx});
    }
}