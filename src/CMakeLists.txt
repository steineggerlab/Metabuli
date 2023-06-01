include_directories(commons)
include_directories(workflow)
include_directories(../lib/prodigal)

add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(metabuli
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        metabuli.cpp
        LocalCommandDeclarations.h util/filter_by_genus.cpp)

mmseqs_setup_derived_target(metabuli)
target_link_libraries(metabuli version)
target_link_libraries(metabuli prodigal)
add_dependencies(metabuli local-generated)

install(TARGETS metabuli DESTINATION bin)
