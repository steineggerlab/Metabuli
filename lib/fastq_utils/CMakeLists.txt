add_library(fastq_utils
        fastq_info.cpp
        hash.cpp
        hash.h
        fastq.cpp
        fastq.h
)

set_target_properties(fastq_utils PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")