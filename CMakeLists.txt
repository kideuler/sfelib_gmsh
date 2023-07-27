cmake_minimum_required(VERSION 3.22)

if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else()
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose build type")
endif()


project(sfelib_gmsh CXX C)
add_library(sfelib_gmsh)
add_subdirectory(include)

set(CMAKE_PREFIX_PATH sfelib_gmsh)

set(CMAKE_CXX_STANDARD 11)
find_library(GMSH_LIB NAMES gmsh)
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()

include_directories(${PROJECT_SOURCE_DIR}/gmsh_install/include)
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(main main.cpp)
target_link_libraries(main PUBLIC ${sfelib_gmsh})
target_link_libraries(main PUBLIC ${GMSH_LIB})

