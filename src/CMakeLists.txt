include_directories(../lib2)
#link_directories(../build/lib2)

add_executable(ADclassifier2
        main.cpp
        BitManipulateMacros.h
        IndexCreator.cpp
        IndexCreator.h
        Kmer.h
        KmerExtractor.cpp
        KmerExtractor.h
        DiffIdxMerger.h
        DiffIdxMerger.cpp
        Mmap.h
        printBinary.cpp
        printBinary.h
        Searcher.cpp
        Searcher.h
        createTargetDB.cpp
        createTargetDB.h common.h fileUtil.cpp fileUtil.h )

target_link_libraries(ADclassifier2 lib2 )