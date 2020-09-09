cmake_minimum_required(VERSION 3.16)
project(ADclassifier2)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(src)
add_subdirectory(lib2)

include_directories(lib/mmseqs)
