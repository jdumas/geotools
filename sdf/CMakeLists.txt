################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)
project(sdf)

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
		mesh_AABB.cpp
		mesh_AABB.h
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

# OpenCL
find_package(OpenCL REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenCL::OpenCL)

# Geogram library
include(geogram)
target_link_libraries(${PROJECT_NAME} PUBLIC geogram)

# Boost.Compute
find_package(Boost 1.61 REQUIRED)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})

# OpenMP
find_package(OpenMP REQUIRED)
target_compile_options(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS})

################################################################################

set(KERNEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/kernels/")
target_compile_definitions(${PROJECT_NAME} PUBLIC -DKERNEL_DIR=\"${KERNEL_DIR}\")
