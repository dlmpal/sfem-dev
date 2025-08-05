#include "fe_factory.hpp"

namespace sfem::fem
{
    std::shared_ptr<NodalFiniteElement>
    create_nodal_element(mesh::CellType cell_type, int order)
    {
        NodalFiniteElement *element = nullptr;

        switch (cell_type)
        {
        case mesh::CellType::point:
            element = new fixed_order::Point;
            break;
        case mesh::CellType::line:
            switch (order)
            {
            case 1:
                element = new fixed_order::Line2;
                break;
            case 2:
                element = new fixed_order::Line3;
                break;
            case 3:
                element = new fixed_order::Line4;
                break;
            default:
                break;
            }
            break;
        case mesh::CellType::triangle:
            switch (order)
            {
            case 1:
                element = new fixed_order::Tri3;
                break;
            case 2:
                element = new fixed_order::Tri6;
                break;
            case 3:
                element = new fixed_order::Tri10;
                break;
            default:
                break;
            }
            break;
        case mesh::CellType::quadrilateral:
            switch (order)
            {
            case 1:
                element = new fixed_order::Quad4;
                break;
            case 2:
                element = new fixed_order::Quad9;
                break;
            case 3:
                element = new fixed_order::Quad16;
                break;
            default:
                break;
            }
            break;
        case mesh::CellType::tetrahedron:
            switch (order)
            {
            case 1:
                element = new fixed_order::Tet4;
                break;
            case 2:
                element = new fixed_order::Tet10;
                break;
            default:
                break;
            }
            break;
        case mesh::CellType::hexahedron:
            switch (order)
            {
            case 1:
                element = new fixed_order::Hex8;
                break;
            default:
                break;
            }
            break;
        default:
            break;
        }

        return std::shared_ptr<NodalFiniteElement>(element);
    }
}