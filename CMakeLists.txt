cmake_minimum_required(VERSION 3.10.2)
project(ADclassifier2)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

add_subdirectory(src)
add_subdirectory(lib2)

include_directories(lib/mmseqs)
