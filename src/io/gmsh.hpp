#pragma once

#include "../mesh/mesh.hpp"
#include <memory>
#include <filesystem>

/// @brief Gmsh file I/O
namespace sfem::io::gmsh
{
    /// @brief Read a Gmsh file
    /// @note Should be in ASCII 2 format
    /// @note Not MPI collective
    std::shared_ptr<mesh::Mesh> read(const std::filesystem::path &filename);
}