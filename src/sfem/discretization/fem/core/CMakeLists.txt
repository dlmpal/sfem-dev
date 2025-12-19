#==============================================================================
add_subdirectory(quadrature)
#==============================================================================
add_subdirectory(elements)
#==============================================================================
add_subdirectory(utils)
#==============================================================================
target_sources(sfem PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/fe_space.cpp
${CMAKE_CURRENT_SOURCE_DIR}/cg_space.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fe_field.cpp
${CMAKE_CURRENT_SOURCE_DIR}/dirichlet_bc.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fe_equation.cpp)