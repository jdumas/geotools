################################################################################
cmake_minimum_required(VERSION 3.3)
project(geotools)
################################################################################

add_subdirectory(convert_mesh)
add_subdirectory(normalize_mesh)
add_subdirectory(poisson_disk)
add_subdirectory(voxmesh)
if(NOT GEOTOOLS_BUILD_ON_TRAVIS)
	add_subdirectory(sdf)
endif()
