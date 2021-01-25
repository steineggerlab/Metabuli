cmake_minimum_required(VERSION 3.10.2)
project(adclassifier2 CXX)
set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake")

set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)
add_subdirectory(lib/mmseqs)
add_subdirectory(src)
#include_directories(lib/prodigal)
include_directories(lib)

add_subdirectory(lib/prodigal)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp main.cpp -o main -lomp")



#include_directories(lib/mmseq)
