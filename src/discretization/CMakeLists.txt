#==============================================================================
add_subdirectory(fem)
#==============================================================================
add_subdirectory(fvm)
#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/function.cpp
${CMAKE_CURRENT_SOURCE_DIR}/coefficient.cpp)