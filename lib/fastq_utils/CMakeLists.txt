add_library(fastq_utils
        fastq_info.cpp
        hash.cpp
        hash.h
        fastq.cpp
        fastq.h
)

find_package(ZLIB REQUIRED)
target_include_directories(fastq_utils PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(fastq_utils PUBLIC ${ZLIB_LIBRARIES})

set_target_properties(fastq_utils PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")
