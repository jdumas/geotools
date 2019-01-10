cmake_minimum_required(VERSION 3.3)
get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
include(geotools)

################################################################################

geotools_import(geogram opencl openmp compute)
geotools_add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME}
	geogram::geogram
	Boost::compute
	OpenCL::OpenCL
	OpenMP::OpenMP_CXX
)

set(KERNEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/kernels/")
target_compile_definitions(${PROJECT_NAME} PUBLIC -DKERNEL_DIR=\"${KERNEL_DIR}\")
