cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(tracer)
set (CMAKE_CXX_STANDARD 14)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake")
include(one-ring-common)

set(HEADER_FILES include/Tracer.h
        include/tracer_c_interface.h)

add_library(tracer SHARED include/Tracer.cpp include/tracer_c_interface.cpp include/TracerImpl.h)
add_library(tracer_static  STATIC include/Tracer.cpp include/tracer_c_interface.cpp include/TracerImpl.h)
if(WIN32)
  set_target_properties(tracer_static PROPERTIES OUTPUT_NAME tracer_static)
else()
  set_target_properties(tracer_static PROPERTIES OUTPUT_NAME tracer)
endif()
add_library(tracer::tracer ALIAS tracer)
add_library(tracer::tracer_static ALIAS tracer_static)
target_include_directories(tracer PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )

target_include_directories(tracer_static PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        )

if(RING_BUILD_TESTING)
    add_subdirectory(tests)
endif()

add_subdirectory(utils)

install(TARGETS tracer tracer_static
        EXPORT tracer
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

install(FILES ${HEADER_FILES} include/tracer_fortran_interface.inc DESTINATION include)

install(EXPORT tracer
        NAMESPACE tracer:: DESTINATION cmake)
install(FILES tracerConfig.cmake
        DESTINATION cmake)


