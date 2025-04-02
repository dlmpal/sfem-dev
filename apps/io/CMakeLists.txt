#=============================================================================
# gmsh-to-sfem
add_executable(gmsh-to-sfem gmsh_to_sfem.cpp)
target_link_libraries(gmsh-to-sfem PRIVATE sfem)
#=============================================================================
# sfem-to-vtk
add_executable(sfem-to-vtk sfem_to_vtk.cpp)
target_link_libraries(sfem-to-vtk PRIVATE sfem)
#==============================================================================
# Installation
include(GNUInstallDirs)
install( 
  TARGETS gmsh-to-sfem
          sfem-to-vtk
  EXPORT sfemTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})