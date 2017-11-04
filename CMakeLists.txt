cmake_minimum_required(VERSION 3.8)
project(Hooke_Jeeves)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp
        include/HJWrapper.h src/HJWrapper.cpp)
add_executable(Hooke_Jeeves ${SOURCE_FILES})
