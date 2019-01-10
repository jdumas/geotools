cmake_minimum_required(VERSION 3.3)
get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
include(geotools)

################################################################################

geotools_import(geogram)
geotools_add_executable(${PROJECT_NAME} main.cpp random.cpp)
target_link_libraries(${PROJECT_NAME} geogram::geogram)
