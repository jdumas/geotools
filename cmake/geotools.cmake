################################################################################

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	get_directory_property(HAS_PARENT PARENT_DIRECTORY)
	if(HAS_PARENT)
		set(CMAKE_BUILD_TYPE "Release" PARENT_SCOPE)
	else()
		set(CMAKE_BUILD_TYPE "Release")
	endif()
endif()

set(GEOTOOLS_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty")
include(GeotoolsDownloadExternal)

# Color output
include(UseColors)

################################################################################

# Eigen
function(geotools_import_eigen)
	if(NOT TARGET Eigen3::Eigen)
		geotools_download_eigen()
		add_library(eigen_eigen INTERFACE)
		add_library(Eigen3::Eigen ALIAS eigen_eigen)
		target_include_directories(eigen_eigen SYSTEM INTERFACE ${GEOTOOLS_EXTERNAL}/eigen)
	endif()
endfunction()

# libigl
function(geotools_import_libigl)
	if(NOT TARGET igl::core)
		geotools_download_libigl()
		list(APPEND CMAKE_MODULE_PATH "${GEOTOOLS_EXTERNAL}/libigl/cmake")
		include(libigl)
	endif()
endfunction()

# Geogram
function(geotools_import_geogram)
	if(NOT TARGET geogram::geogram)
		geotools_download_geogram()
		set(GEOGRAM_SEARCH_PATHS "${GEOTOOLS_EXTERNAL}/geogram")
		include(geogram)
	endif()
endfunction()

# OpenCL
function(geotools_import_opencl)
	if(NOT TARGET OpenCL::OpenCL)
		cmake_minimum_required(VERSION 3.7)
		find_package(OpenCL REQUIRED)
	endif()
endfunction()

# Boost.Compute
function(geotools_import_compute)
	# if(NOT TARGET Boost::compute)
	# 	find_package(Boost 1.61 REQUIRED COMPONENTS compute QUIET)
	# 	if(NOT TARGET Boost::compute)
	# 		# When CMake and Boost versions are not in sync, imported targets may not be available... (sigh)
	# 		add_library(boost_compute INTERFACE)
	# 		add_library(Boost::compute ALIAS boost_compute)
	# 		target_include_directories(boost_compute SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
	# 			target_link_libraries(boost_compute INTERFACE ${Boost_LIBRARIES})
	# 	endif()
	# endif()
	if(NOT TARGET Boost::compute)
		geotools_download_compute()
		add_library(boost_compute INTERFACE)
		add_library(Boost::compute ALIAS boost_compute)
		target_include_directories(boost_compute SYSTEM INTERFACE "${GEOTOOLS_EXTERNAL}/compute/include")
	endif()
endfunction()

# OpenMP
function(geotools_import_openmp)
	if(NOT TARGET OpenMP::OpenMP_CXX)
		cmake_minimum_required(VERSION 3.10)
		find_package(OpenMP REQUIRED)
	endif()
endfunction()

################################################################################

# Add executable
function(geotools_add_executable name)
	add_executable(${name} ${ARGN})
	target_link_libraries(${name} colors::colors)

	# Use C++11
	target_compile_features(${name} PUBLIC cxx_std_11)

	# Output folder
	set_target_properties(${name}
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endfunction()

# All
set(__geotools_import_dir "${CMAKE_CURRENT_LIST_DIR}")
function(geotools_import)
	foreach(NAME IN ITEMS ${ARGN})
		set(__import_file "${CMAKE_CURRENT_BINARY_DIR}/geotools_import_${NAME}.cmake")
		configure_file("${__geotools_import_dir}/GeotoolsImport.cmake.in" "${__import_file}" @ONLY)
		include("${__import_file}")
	endforeach()

	if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		foreach(config ${CMAKE_CONFIGURATION_TYPES})
			string(TOUPPER ${config} config)
			string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
			string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
		endforeach()
	endif()
endfunction()
