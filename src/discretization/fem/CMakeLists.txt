#==============================================================================
add_subdirectory(quadrature)
#==============================================================================
add_subdirectory(elements)
#==============================================================================
add_subdirectory(utils)
#==============================================================================
add_subdirectory(kernels)
#==============================================================================
target_sources(sfem PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/fe_space.cpp
${CMAKE_CURRENT_SOURCE_DIR}/cg_space.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fe_function.cpp
${CMAKE_CURRENT_SOURCE_DIR}/dirichlet_bc.cpp)