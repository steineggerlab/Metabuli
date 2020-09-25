include_directories(commons)
add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)

include_directories(../lib2)
include_directories(../lib)
#link_directories(../build/lib2)

add_executable(ADclassifier2
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        main.cpp
        )

target_link_libraries(ADclassifier2 version)
