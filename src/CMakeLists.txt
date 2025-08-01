include_directories(commons)
include_directories(workflow)
include_directories(util)
include_directories(benchmark)
include_directories(../lib/prodigal)
include_directories(../lib/fasta_validator)
include_directories(../lib/fastq_utils)

add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)
add_subdirectory(benchmark)

add_executable(metabuli
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        ${benchmark_source_files}
        metabuli.cpp
        LocalCommandDeclarations.h 
        util/filter_by_genus.cpp
        MetabuliBase.cpp)

mmseqs_setup_derived_target(metabuli)
target_link_libraries(metabuli version)
target_link_libraries(metabuli prodigal)
target_link_libraries(metabuli fasta_validator)
target_link_libraries(metabuli fastq_utils)
add_dependencies(metabuli local-generated)

install(TARGETS metabuli DESTINATION bin)
