#==============================================================================
add_subdirectory(utils)
#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/cell.cpp
${CMAKE_CURRENT_SOURCE_DIR}/region.cpp
${CMAKE_CURRENT_SOURCE_DIR}/topology.cpp
${CMAKE_CURRENT_SOURCE_DIR}/partition.cpp
${CMAKE_CURRENT_SOURCE_DIR}/mesh.cpp)