#=============================================================================
# simple
add_executable(simple simple.cpp)
target_link_libraries(simple PRIVATE sfem)
#==============================================================================
# Installation
include(GNUInstallDirs)
install( 
  TARGETS simple
  EXPORT sfemTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})