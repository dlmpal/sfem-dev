#pragma once

#include <sfem/geo/vec3.hpp>
#include <sfem/base/error.hpp>
#include <span>

#define SFEM_BAD_CELL_ERROR(cell_type)                                                    \
    {                                                                                     \
        auto msg = std::format("Cell type {} is invalid\n", static_cast<int>(cell_type)); \
        SFEM_ERROR(msg);                                                                  \
    }

namespace sfem::mesh
{
    /// @brief Available cell types
    enum class CellType : int
    {
        point = 0,
        line = 1,
        triangle = 2,
        quadrilateral = 3,
        tetrahedron = 4,
        hexahedron = 5,
        prism = 6,
        n_cell_types = 7
    };

    /// @brief Get the string representation for given cell type
    std::string cell_type_str(CellType cell_type);

    /// @brief Get the dimension for a given cell type
    int cell_dim(CellType cell_type);

    /// @brief Get the number of nodes for a given cell type
    int cell_num_nodes(CellType cell_type);

    /// @brief Get the number of edges for a given cell type
    int cell_num_edges(CellType cell_type);

    /// @brief Get the number of faces for a given cell type
    int cell_num_faces(CellType cell_type);

    /// @brief Get the ordering of the edge nodes for a specific cell type
    std::array<int, 2> cell_edge_ordering(CellType cell_type, int edge_idx);

    /// @brief Get the face CellType for a specific cell type and face index
    CellType cell_face_type(CellType cell_type, int face_idx);

    /// @brief Get the ordering of the face nodes for a specific cell type
    std::array<int, 4> cell_face_ordering(CellType cell_type, int face_idx);

    /// @brief Simple datatype that stores the region tag
    /// and type for a cell
    struct Cell
    {
        int tag = -1;
        CellType type = CellType::point;
    };
}