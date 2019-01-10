################################################################################
cmake_minimum_required(VERSION 3.3)
project(geotools)
################################################################################

add_subdirectory(convert_mesh)
add_subdirectory(normalize_mesh)
add_subdirectory(poisson_disk)
add_subdirectory(voxmesh)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	add_subdirectory(sdf)
endif()
