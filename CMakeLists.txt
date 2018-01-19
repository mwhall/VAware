# configures building lib
cmake_minimum_required(VERSION 2.6)
project(VAware)

# versioning
set(VAware_VERSION_MAJOR 0)
set(VAware_VERSION_MINOR 0)
set(VAware_VERSION_PATCH 1)
set(VAware_VERSION "${VAware_VERSION_MAJOR}.${VAware_VERSION_MINOR}.${VAware_VERSION_PATCH}")

# config version headers
configure_file (
    "${PROJECT_SOURCE_DIR}/vaware_config.h.in"
    "${PROJECT_BINARY_DIR}/vaware_config.h"
    )
include_directories("${PROJECT_BINARY_DIR}")

# find dependencies
find_package(Threads REQUIRED)
find_package(SeqAn REQUIRED)

# build
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
add_executable(vaware vaware.cpp)
target_link_libraries(vaware ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
