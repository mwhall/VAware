# configures building lib
cmake_minimum_required(VERSION 2.6)
find_package(Threads)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
add_executable(vaware vaware.cpp)
target_link_libraries(vaware ${SEQAN_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
