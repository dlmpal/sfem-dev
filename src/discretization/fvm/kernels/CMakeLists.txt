#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/laplacian.cpp
${CMAKE_CURRENT_SOURCE_DIR}/convection.cpp
${CMAKE_CURRENT_SOURCE_DIR}/transient.cpp
${CMAKE_CURRENT_SOURCE_DIR}/source.cpp)