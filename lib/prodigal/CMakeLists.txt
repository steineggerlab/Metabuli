add_library(prodigal
        bitmap.cpp
        bitmap.h
        dprog.cpp
        dprog.h
        fptr.h
        gene.cpp
        gene.h
        metagenomic.cpp
        metagenomic.h
        node.cpp
        node.h
        prodigalsequence.cpp
        prodigalsequence.h
        training.cpp
        training.h)

set_target_properties(prodigal PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")
