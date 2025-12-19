#include <sfem/sfem.hpp>

using namespace sfem;

int main(int argc, char **argv)
{
    initialize(argc, argv, "sfem-to-vtk");
    std::string mesh_dir = argv[1];
    std::string vtk_filename = argv[2];
    auto mesh = io::read_mesh(mesh_dir);
    io::vtk::write(vtk_filename, *mesh);
    return 0;
}