#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/flux_function.cpp
${CMAKE_CURRENT_SOURCE_DIR}/numerical_flux.cpp
${CMAKE_CURRENT_SOURCE_DIR}/euler.cpp)