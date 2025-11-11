cmake_minimum_required(VERSION 3.16)
project(conv-diff)

set(CMAKE_CXX_COMPILER $ENV{MPICXX})
set(sfem_DIR $ENV{SFEM_DEV_INSTALL_DIR}/lib/cmake/sfem)
find_package(sfem REQUIRED)

add_executable(conv-diff ../main.cpp)
target_link_libraries(conv-diff sfem::sfem)