cmake_minimum_required(VERSION 3.22)

project(SFE_GMSH CXX C)
add_subdirectory(src)

set(CMAKE_PREFIX_PATH SFE_GMSH)

set(CMAKE_CXX_STANDARD 11)
find_library(GMSH_LIB NAMES gmsh PATHS "${PROJECT_SOURCE_DIR}/gmsh_install/lib")
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

find_path(GMSH_INC gmsh.h PATHS "${PROJECT_SOURCE_DIR}/gmsh_install/include")
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()



include_directories(${PROJECT_SOURCE_DIR}/gmsh_install/include)
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE SFE_GMSH)
target_link_libraries(main PRIVATE ${GMSH_LIB})

