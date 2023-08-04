set(HEADER_LIST "${PROJECT_SOURCE_DIR}/include/sfelib_gmsh.hpp")

add_library(SFE_GMSH sfelib_gmsh_kernals.cpp ${HEADER_LIST})


target_include_directories(SFE_GMSH PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SFE_GMSH PUBLIC ${PROJECT_SOURCE_DIR}/gmsh_install/include)



source_group(
TREE "${PROJECT_SOURCE_DIR}/include"

PREFIX "Header Files"

FILES ${HEADER_LIST})