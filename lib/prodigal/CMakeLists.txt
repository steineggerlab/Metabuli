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

mmseqs_setup_derived_target(prodigal)