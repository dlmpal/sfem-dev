#=============================================================================
# cart-mesh
set(argparse_DIR ${ARGPARSE_DIR}/lib/cmake/argparse)
find_package(argparse REQUIRED)

add_executable(cart-mesh cart-mesh.cpp)
target_link_libraries(cart-mesh PRIVATE sfem argparse::argparse)
#==============================================================================
# Installation
include(GNUInstallDirs)
install( 
  TARGETS cart-mesh
  EXPORT sfemTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})