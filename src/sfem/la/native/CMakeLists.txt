#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/dense_matrix.cpp
${CMAKE_CURRENT_SOURCE_DIR}/dense_matrix_utils.cpp
${CMAKE_CURRENT_SOURCE_DIR}/sparsity.cpp
${CMAKE_CURRENT_SOURCE_DIR}/vector.cpp
${CMAKE_CURRENT_SOURCE_DIR}/sparse_matrix.cpp
${CMAKE_CURRENT_SOURCE_DIR}/setval_utils.cpp
${CMAKE_CURRENT_SOURCE_DIR}/linear_system.cpp)
#==============================================================================
add_subdirectory(linear_solvers)