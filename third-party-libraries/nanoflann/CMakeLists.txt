cmake_minimum_required(VERSION 3.6 FATAL_ERROR)
project(nanoflann LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
  nanoflann_upstream
  GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
  GIT_TAG        v1.3.0
  GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(nanoflann_upstream)
if(NOT nanoflann_upstream_POPULATED)
  FetchContent_Populate(nanoflann_upstream)
endif()

add_library(nanoflann INTERFACE)
add_library(nanoflann::nanoflann ALIAS nanoflann)
target_include_directories(nanoflann INTERFACE
  $<BUILD_INTERFACE:${nanoflann_upstream_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

install(TARGETS nanoflann EXPORT nanoflann)
install(EXPORT nanoflann
  NAMESPACE nanoflann::
  DESTINATION cmake)
install(FILES nanoflannConfig.cmake DESTINATION cmake)
