cmake_minimum_required(VERSION 3.16)
project(elasticity-3d)

set(CMAKE_CXX_COMPILER $ENV{MPICXX})

set(sfem_DIR $ENV{SFEM_DEV_INSTALL_DIR}/lib/cmake/sfem)
find_package(sfem REQUIRED)

set(argparse_DIR /home/user/local/argparse/lib/cmake/argparse)
find_package(argparse REQUIRED)

add_executable(elasticity-3d ../main.cpp)
target_link_libraries(elasticity-3d sfem::sfem argparse::argparse)