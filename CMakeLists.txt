cmake_minimum_required(VERSION 3.8)
project(Hooke_Jeeves)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES src/main.cpp
        include/HJWrapper.h src/HJWrapper.cpp)

add_subdirectory(lib/googletest)
include_directories(lib/googletest/googletest/include)
include_directories(lib/googletest/googlemock/include)

add_executable(Hooke_Jeeves_exe ${SOURCE_FILES})
add_library(Hooke_Jeeves_lib SHARED ${SOURCE_FILES})
target_link_libraries(Hooke_Jeeves_exe gtest gtest_main)

add_subdirectory(test)
