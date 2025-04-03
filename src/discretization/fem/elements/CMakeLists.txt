#==============================================================================
target_sources(sfem PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/fe.cpp
${CMAKE_CURRENT_SOURCE_DIR}/nodal_fe.cpp
${CMAKE_CURRENT_SOURCE_DIR}/point.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fixed_order_line.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fixed_order_triangle.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fixed_order_quad.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fixed_order_tet.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fixed_order_hex.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fe_factory.cpp)