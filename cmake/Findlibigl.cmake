# Find libigl(http://igl.ethz.ch/projects/libigl/)
# The following variables are set
#
# LIBIGL_FOUND
# LIBIGL_INCLUDE_DIRS
#
# It searches the environment variable $LIBIGL_PATH

if(PREFER_LOCAL_LIBS)
	find_path(LIBIGL_INCLUDE
			igl/igl_inline.h
			PATHS ${THIRD_PARTY_DIR}/libigl_
			PATH_SUFFIXES include
			NO_DEFAULT_PATH
	)
endif()

# If nothing is found, search against but include system paths
find_path(LIBIGL_INCLUDE
		igl/igl_inline.h
		HINTS
			ENV LIBIGL_PATH
		PATHS
			${THIRD_PARTY_DIR}/libigl
			"C:/Program Files/libigl/"
		PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBIGL DEFAULT_MSG LIBIGL_INCLUDE)

if(LIBIGL_FOUND)
	set(LIBIGL_INCLUDE_DIRS ${LIBIGL_INCLUDE})
endif()

mark_as_advanced(LIBIGL_INCLUDE_DIRS)

################################################################################

# libigl options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
# option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)

# Adding libigl: choose the path to your local copy libigl
# This is going to compile everything you requested
#message(FATAL_ERROR "${PROJECT_SOURCE_DIR}/../libigl/cmake")
if(NOT TARGET igl)
	set(LIBIGL_SOURCE_DIR ${LIBIGL_INCLUDE_DIRS})

	# Target igl
	if(LIBIGL_USE_STATIC_LIBRARY)
		file(GLOB SOURCES_IGL
			"${LIBIGL_SOURCE_DIR}/igl/*.cpp"
			"${LIBIGL_SOURCE_DIR}/igl/copyleft/*.cpp")
		add_library(igl STATIC ${SOURCES_IGL})
		target_include_directories(igl SYSTEM PUBLIC ${LIBIGL_SOURCE_DIR})
		target_compile_definitions(igl PUBLIC -DIGL_STATIC_LIBRARY)

		# Use C++11
		set_target_properties(igl PROPERTIES CXX_STANDARD 11)
		set_target_properties(igl PROPERTIES CXX_STANDARD_REQUIRED ON)

		# Eigen3 library
		find_package(Eigen3 REQUIRED)
		target_include_directories(igl SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})

		# Check program behavior
		find_package(Sanitizers)
		add_sanitizers(igl)
	else()
		# Interface target for libigl
		add_library(igl INTERFACE)
		target_include_directories(igl SYSTEM INTERFACE ${LIBIGL_SOURCE_DIR})

		# Eigen3 library
		find_package(Eigen3 REQUIRED)
		target_include_directories(igl SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
	endif()
endif()
