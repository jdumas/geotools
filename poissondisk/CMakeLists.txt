# Project
CMake_Minimum_Required (VERSION 3.0)
Project(poisson_disk)

# Default build mode
If (NOT CMAKE_BUILD_TYPE)
	Message (STATUS "No build type selected, default to Release")
	Set (CMAKE_BUILD_TYPE "Release")
Endif ()

# Executable
Add_Executable(poisson_disk main.cpp random.cpp)

# Use C++11
Set_Target_Properties(poisson_disk PROPERTIES CXX_STANDARD 11)
Set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check program behavior in Debug mode
If (CMAKE_BUILD_TYPE MATCHES Debug)
	Target_Compile_Options (poisson_disk INTERFACE "-fsanitize=address")
Endif ()
