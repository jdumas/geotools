################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)
project(voxmesh)

################################################################################

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

################################################################################

# Main files
add_executable(${PROJECT_NAME}
		main.cpp
		octree.cpp

		attributes.h
		common.h
		octree.h
)

# C++11
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

# Output binary location
set_target_properties(${PROJECT_NAME}
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

################################################################################
# 3rd party libraries
################################################################################

# Eigen library
find_package(Eigen3 REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

# Geogram library
include(geogram)
target_link_libraries(${PROJECT_NAME} PUBLIC geogram)
