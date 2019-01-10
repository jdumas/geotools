################################################################################

include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
set(GEOTOOLS_EXTRA_OPTIONS TLS_VERIFY OFF)
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
	list(APPEND GEOTOOLS_EXTRA_OPTIONS GIT_CONFIG advice.detachedHead=false)
endif()

option(GEOTOOLS_SKIP_DOWNLOAD "Skip downloading external libraries" OFF)

# Shortcut functions
function(geotools_download_project name)
	if(NOT GEOTOOLS_SKIP_DOWNLOAD)
		download_project(
			PROJ         ${name}
			SOURCE_DIR   "${GEOTOOLS_EXTERNAL}/${name}"
			DOWNLOAD_DIR "${GEOTOOLS_EXTERNAL}/.cache/${name}"
			QUIET
			${GEOTOOLS_EXTRA_OPTIONS}
			${ARGN}
		)
	endif()
endfunction()

################################################################################

## Eigen
function(geotools_download_eigen)
	geotools_download_project(eigen
		URL           http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
		URL_MD5       f2a417d083fe8ca4b8ed2bc613d20f07
	)
endfunction()

## libigl
function(geotools_download_libigl)
	geotools_download_project(libigl
		GIT_REPOSITORY https://github.com/libigl/libigl.git
		GIT_TAG        7895811be45816ba399365a13487305660f152a7
	)
endfunction()

## geogram
function(geotools_download_geogram)
	geotools_download_project(geogram
		GIT_REPOSITORY https://github.com/alicevision/geogram.git
		GIT_TAG        v1.6.9
	)
endfunction()

# Boost.Compute
function(geotools_download_compute)
	geotools_download_project(compute
		GIT_REPOSITORY https://github.com/boostorg/compute
		GIT_TAG        9189a761b79fcd4be2f38158b9cad164bac22fa2
	)
endfunction()
