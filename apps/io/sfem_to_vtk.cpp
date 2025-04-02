#include "sfem.hpp"

int main(int argc, char **argv)
{
    sfem::Application::instance(argc, argv, "sfem-to-vtk");
    std::string mesh_dir = argv[1];
    std::string vtk_filename = argv[2];
    auto mesh = sfem::io::read_mesh(mesh_dir);
    sfem::io::vtk::write(vtk_filename, *mesh);
    return 0;
}