cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
# keep old policy for setting implicit link libraries
# zlib causes issues in static builds otherwise
cmake_policy(SET CMP0060 OLD)

project(metabuli LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 14)

# Module paths
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake"
)

# Default to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# MMseqs framework setup
set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)

# Third-party: libdeflate (with CMake)
add_subdirectory(lib/libdeflate)

# Third-party: libisal (no CMakeLists, build manually)
# Collect all C/C++ sources, exclude architecture-specific files
file(GLOB_RECURSE ISAL_SRC_TMP
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/libisal/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/libisal/*.cc"
)
# Filter out architecture-specific sources (e.g., aarch64, ppc64le)
set(ISAL_SRC "")
foreach(_src IN LISTS ISAL_SRC_TMP)
    if(NOT _src MATCHES ".*(aarch64|ppc64le|vsx|riscv|tests|examples).*" )
        list(APPEND ISAL_SRC ${_src})
    endif()
endforeach()
if (NOT ISAL_SRC)
    message(FATAL_ERROR "No libisal sources found (after filtering) in lib/libisal/src")
endif()
add_library(isal STATIC ${ISAL_SRC})
target_include_directories(isal PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/libisal/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/libisal" 
)




# fastp: no CMakeLists, build manually as static library
file(GLOB_RECURSE FASTP_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/fastp/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/fastp/*.c"
)
if (NOT FASTP_SRC)
    message(FATAL_ERROR "No fastp sources found in lib/fastp/src")
endif()
add_library(fastp STATIC ${FASTP_SRC})
target_include_directories(fastp PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/fastp"
)

# Link deflate and isal to fastp
target_link_libraries(fastp PUBLIC deflate isal)


add_subdirectory(lib/mmseqs)
add_subdirectory(lib/prodigal)

add_subdirectory(src)
add_subdirectory(data)

include_directories(lib)
# # Include top-level lib headers for all targets
# include_directories(
#     "${CMAKE_CURRENT_SOURCE_DIR}/lib/libdeflate/include"
#     "${CMAKE_CURRENT_SOURCE_DIR}/lib/libisal/include"
#     "${CMAKE_CURRENT_SOURCE_DIR}/lib/libfastp/include"
#     "${CMAKE_CURRENT_SOURCE_DIR}/lib"
# )



