#==============================================================================
target_sources(sfem PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/petsc_vec.cpp
${CMAKE_CURRENT_SOURCE_DIR}/petsc_mat.cpp
${CMAKE_CURRENT_SOURCE_DIR}/petsc_ksp.cpp
${CMAKE_CURRENT_SOURCE_DIR}/petsc.cpp
${CMAKE_CURRENT_SOURCE_DIR}/petsc_linear_system.cpp)