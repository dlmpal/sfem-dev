#==============================================================================
add_subdirectory(solid_mechanics)
#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/mass.cpp
${CMAKE_CURRENT_SOURCE_DIR}/diffusion.cpp)