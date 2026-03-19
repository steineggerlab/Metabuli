include_directories(commons)
include_directories(workflow)
include_directories(util)
include_directories(benchmark)
include_directories(uniref)
include_directories(read-group)
include_directories(../lib/prodigal)
include_directories(../lib/fasta_validator)
include_directories(../lib/fastq_utils)
include_directories(../lib/yxml)
include_directories(../lib/json)

add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)
add_subdirectory(benchmark)
add_subdirectory(uniref)
add_subdirectory(read-group)

add_executable(metabuli
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        ${benchmark_source_files}
        ${uniref_source_files}
        ${read_group_source_files}
        metabuli.cpp
        LocalCommandDeclarations.h 
        util/filter_by_genus.cpp
        MetabuliBase.cpp)

mmseqs_setup_derived_target(metabuli)

if (HAVE_SSE2)
    target_compile_options(metabuli PRIVATE -msse4.1 -mno-avx512f -mno-avx512vl)
elseif (HAVE_AVX2)
    target_compile_options(metabuli PRIVATE -mavx2 -mno-avx512f -mno-avx512vl)
else()
    target_compile_options(metabuli PRIVATE -mno-avx512f -mno-avx512vl)
endif()

target_link_libraries(metabuli version)
target_link_libraries(metabuli prodigal)
target_link_libraries(metabuli fasta_validator)
target_link_libraries(metabuli fastq_utils)
target_link_libraries(metabuli yxml)
target_link_libraries(metabuli jsonlib)
add_dependencies(metabuli local-generated)

install(TARGETS metabuli DESTINATION bin)

# =====================================================================
# THE AVX-512 SLEDGEHAMMER 
# =====================================================================

# 1. Rip -march=native out of the global CMake flags
string(REPLACE "-march=native" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# 2. Explicitly disable AVX-512 and its Vector Length (VL) extensions.
# We must disable 'vl' because that is what triggers avx512vlintrin.h.
target_compile_options(metabuli PRIVATE "-mno-avx512f" "-mno-avx512vl" "-mno-avx512dq" "-mno-avx512cd" "-mno-avx512bw")

# 3. Disarm MMseqs's hidden interface flags (This is what bypassed the last fix)
if (TARGET mmseqs)
    get_target_property(_mmseqs_opts mmseqs INTERFACE_COMPILE_OPTIONS)
    if (_mmseqs_opts)
        list(REMOVE_ITEM _mmseqs_opts "-march=native")
        set_target_properties(mmseqs PROPERTIES INTERFACE_COMPILE_OPTIONS "${_mmseqs_opts}")
    endif()
endif()

if (TARGET mmseqs_core)
    get_target_property(_mmseqs_core_opts mmseqs_core INTERFACE_COMPILE_OPTIONS)
    if (_mmseqs_core_opts)
        list(REMOVE_ITEM _mmseqs_core_opts "-march=native")
        set_target_properties(mmseqs_core PROPERTIES INTERFACE_COMPILE_OPTIONS "${_mmseqs_core_opts}")
    endif()
endif()
# =====================================================================