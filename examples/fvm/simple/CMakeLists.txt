cmake_minimum_required(VERSION 3.16)
project(simple)

set(CMAKE_CXX_COMPILER $ENV{MPICXX})
set(sfem_DIR $ENV{SFEM_DEV_INSTALL_DIR}/lib/cmake/sfem)
set(argparse_DIR $ENV{ARGPARSE_DIR}/lib/cmake/argparse)
find_package(sfem REQUIRED)
find_package(argparse REQUIRED)

add_executable(simple ../main.cpp)
target_link_libraries(simple sfem::sfem argparse::argparse)