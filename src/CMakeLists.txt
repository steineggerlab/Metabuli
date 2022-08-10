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
        LocalCommandDeclarations.h workflow/inclusion_test.cpp workflow/exclusion_test.cpp workflow/inclusion_test_hiv.cpp workflow/exclusion_test_hiv.cpp commons/Query.cpp commons/Query.h commons/Genus.cpp commons/Genus.h commons/Species.cpp commons/Species.h commons/Match.h commons/common.cpp)

mmseqs_setup_derived_target(metabuli)
target_link_libraries(metabuli version)
target_link_libraries(metabuli prodigal)



install(TARGETS metabuli DESTINATION bin)
