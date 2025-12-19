#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/constitutive.cpp
${CMAKE_CURRENT_SOURCE_DIR}/strain.cpp
${CMAKE_CURRENT_SOURCE_DIR}/stress.cpp
${CMAKE_CURRENT_SOURCE_DIR}/linear_elasticity.cpp
${CMAKE_CURRENT_SOURCE_DIR}/pressure_load.cpp)