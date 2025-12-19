#==============================================================================
target_sources(sfem PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
${CMAKE_CURRENT_SOURCE_DIR}/legacy.cpp
${CMAKE_CURRENT_SOURCE_DIR}/xml_utils.cpp
${CMAKE_CURRENT_SOURCE_DIR}/xml.cpp
${CMAKE_CURRENT_SOURCE_DIR}/vtk.cpp)