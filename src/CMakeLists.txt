include_directories(commons)
include_directories(workflow)
include_directories(../lib/prodigal)
include_directories(../lib/csv-parser)
include_directories(../lib/string-view-lite)

add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(metabuli
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        metabuli.cpp
        LocalCommandDeclarations.h)

mmseqs_setup_derived_target(metabuli)
target_link_libraries(metabuli version)
target_link_libraries(metabuli prodigal)
target_link_libraries(metabuli csv)
target_link_libraries(metabuli string_view_lite)

install(TARGETS metabuli DESTINATION bin)
