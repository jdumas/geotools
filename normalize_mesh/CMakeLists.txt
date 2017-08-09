################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)
project(normalize_mesh)

################################################################################

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

################################################################################

# Main files
add_executable(${PROJECT_NAME} main.cpp)

# Geogram library
include(geogram)
target_link_libraries(${PROJECT_NAME} geogram)

# Target binary folder
set_target_properties(${PROJECT_NAME}
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
