add_library(yxml
    yxml.h
    yxml.c
)

set_target_properties(yxml PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")