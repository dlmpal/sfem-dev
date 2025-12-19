#==============================================================================
# Top level CMakeLists.txt for SFEM
cmake_minimum_required(VERSION 3.16)
project(SFEM VERSION "1.0.0")
#==============================================================================
# Build type
if(NOT DEFINED CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE RELEASE)
endif()

## Relase
if (CMAKE_BUILD_TYPE MATCHES RELEASE)
        set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic -Wconversion -Wno-comment")
endif()

## Debug
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
        set(CMAKE_CXX_FLAGS "-g -Wall -Wextra -Wpedantic -Wconversion")
endif()
#==============================================================================
# Use the build/link RPATH also for installation 
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#==============================================================================
# Create the core library
add_subdirectory(src/sfem)

# Append the library installation location to the RPATH
# This is useful for the installation of downstream targets that link to SFEM 
list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
#==============================================================================
# Create the applications
option(WITH_APPS ON)
if(${WITH_APPS})
        add_subdirectory(apps)
endif()
#==============================================================================
# Installation
install(EXPORT sfemTargets
        FILE sfemTargets.cmake
        NAMESPACE sfem::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sfem)

## Create package config
include(CMakePackageConfigHelpers)
configure_package_config_file( 
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in 
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/sfemConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sfem)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/sfemConfig.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sfem)