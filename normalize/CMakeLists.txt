cmake_minimum_required(VERSION 2.8.11)
project(normalize_mesh)

if (NOT CMAKE_BUILD_TYPE)
	message (STATUS "No build type selected, default to Release")
	set (CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Geogram REQUIRED)

add_executable(normalize_mesh main.cpp)
target_link_libraries(normalize_mesh Geogram::geogram)
