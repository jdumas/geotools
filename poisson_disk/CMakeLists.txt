# Project
cmake_minimum_required(VERSION 3.0)
project(poisson_disk)

# Default build mode
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Executable
add_executable(${PROJECT_NAME} main.cpp random.cpp)

# Use C++11
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

# Check program behavior in Debug mode
if(CMAKE_BUILD_TYPE MATCHES Debug)
	target_compile_options(${PROJECT_NAME} INTERFACE "-fsanitize=address")
endif()

# Output directory for binaries
set_target_properties(${PROJECT_NAME}
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
