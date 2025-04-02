#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/native.cpp
${CMAKE_CURRENT_SOURCE_DIR}/gmsh.cpp
${CMAKE_CURRENT_SOURCE_DIR}/vtk.cpp)