#==============================================================================
add_subdirectory(utils)
#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/fv_space.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fv_function.cpp
${CMAKE_CURRENT_SOURCE_DIR}/diffusion.cpp
${CMAKE_CURRENT_SOURCE_DIR}/hyperbolic.cpp)