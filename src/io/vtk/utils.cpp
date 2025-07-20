#include "utils.hpp"
#include <map>

namespace sfem::io::vtk
{
    //=============================================================================
    int cell_type_to_vtk(mesh::CellType cell_type, int order)
    {
        static const std::map<std::pair<mesh::CellType, int>, int> to_vtk = {
            {{mesh::CellType::point, 1}, 1},

            {{mesh::CellType::line, 1}, 3},
            {{mesh::CellType::line, 2}, 21},
            {{mesh::CellType::line, 3}, 68},

            {{mesh::CellType::triangle, 1}, 5},
            {{mesh::CellType::triangle, 2}, 22},
            {{mesh::CellType::triangle, 3}, 69},

            {{mesh::CellType::quad, 1}, 9},
            {{mesh::CellType::quad, 2}, 23},
            {{mesh::CellType::quad, 3}, 70},

            {{mesh::CellType::tet, 1}, 10},
            {{mesh::CellType::tet, 2}, 24},
            {{mesh::CellType::tet, 3}, 71},

            {{mesh::CellType::hex, 1}, 12},
            {{mesh::CellType::hex, 2}, 25},
            {{mesh::CellType::hex, 3}, 72},

            {{mesh::CellType::prism, 1}, 13}

        };

        if (!to_vtk.contains({cell_type, order}))
        {
            SFEM_ERROR(std::format("Cell type {} of degree {} cannot be converted to VTK equivalent\n",
                                   static_cast<int>(cell_type), order));
        }

        return to_vtk.at({cell_type, order});
    }
    //=============================================================================
    void cell_node_ordering_to_vtk(int vtk_type, std::span<int> nodes)
    {
        // Second-order tetrahedron
        if (vtk_type == 24)
        {
            std::swap(nodes[8], nodes[9]);
        }
    }
}