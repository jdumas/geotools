# CMake
cmake_minimum_required(VERSION 3.0)
project(convert_mesh)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

# Main app
add_executable(${PROJECT_NAME} main.cpp)

# libigl
find_package(libigl REQUIRED)
target_link_libraries(${PROJECT_NAME} igl)

# Geogram
find_package(Geogram REQUIRED)
target_link_libraries(${PROJECT_NAME} Geogram::geogram)

set_target_properties(${PROJECT_NAME}
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
