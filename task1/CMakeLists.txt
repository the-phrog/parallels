cmake_minimum_required(VERSION 3.10)
project(SinArraySum)

option(USE_DOUBLE "Use double instead of float" OFF)

add_executable(sin_array_sum main.cpp)

if(USE_DOUBLE)
    target_compile_definitions(sin_array_sum PRIVATE USE_DOUBLE)
endif()

set_target_properties(sin_array_sum PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
