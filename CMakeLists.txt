cmake_minimum_required(VERSION 3.6)
project(Hooke_Jeeves)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp HookeJeevesWrapper.cpp HookeJeevesWrapper.h)
add_executable(Hooke_Jeeves ${SOURCE_FILES})

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(Hooke_Jeeves PUBLIC Debug=0)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_definitions(Hooke_Jeeves PUBLIC Debug=0)
endif()