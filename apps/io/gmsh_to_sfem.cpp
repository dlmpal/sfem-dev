// Convert a gmsh file to native sfem mesh format
// .msh2 files are suported
// The program is executed as follows:
//   gmshToSfem $gmsh_filename $mesh_directory

#include <sfem/sfem.hpp>

using namespace sfem;

int main(int argc, char **argv)
{
    initialize(argc, argv, "gmsh-to-sfem");
    std::string gmsh_filename = argv[1];
    std::string mesh_dir = argv[2];
    auto mesh = io::gmsh::read(gmsh_filename);
    io::write_mesh(mesh_dir, *mesh);
    return 0;
}