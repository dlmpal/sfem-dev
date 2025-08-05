#include "utils.hpp"
#include <map>

namespace sfem::io::gmsh
{
    //=============================================================================
    mesh::CellType gmsh_type_to_native(int gmsh_type, int gmsh_idx)
    {
        static const std::map<int, mesh::CellType> from_gmsh =
            {
                {15, mesh::CellType::point},
                {1, mesh::CellType::line},
                {2, mesh::CellType::triangle},
                {3, mesh::CellType::quadrilateral},
                {4, mesh::CellType::tetrahedron},
                {5, mesh::CellType::hexahedron},
                {6, mesh::CellType::prism},
            };

        if (from_gmsh.contains(gmsh_type) == false)
        {
            SFEM_ERROR(std::format("Gmsh element {} has unsupported type: {}\n", gmsh_idx, gmsh_type));
        }

        return from_gmsh.at(gmsh_type);
    }
}