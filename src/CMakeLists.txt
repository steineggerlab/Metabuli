include_directories(commons)
include_directories(workflow)
add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(adclassifier
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        adclassifier.cpp
        LocalCommandDeclarations.h
        )
mmseqs_setup_derived_target(adclassifier)
target_link_libraries(adclassifier version)


install(TARGETS adclassifier DESTINATION bin)
