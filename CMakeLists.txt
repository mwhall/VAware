# configures building lib

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
add_executable(vaware vaware.cpp)
target_link_libraries(vaware ${SEQAN_LIBRARIES})
