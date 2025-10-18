#==============================================================================
add_subdirectory(kernels)
#==============================================================================
add_subdirectory(hyperbolic)
#==============================================================================
add_subdirectory(utils)
#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/fv_space.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fv_field.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fv_coefficient.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fv_bc.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fv_gradient.cpp)