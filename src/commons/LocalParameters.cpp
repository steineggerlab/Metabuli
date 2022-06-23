#include "LocalParameters.h"

LocalParameters::LocalParameters() :
        Parameters(),
        PARAM_GTDB_OR_NCBI(PARAM_GTDB_OR_NCBI_ID,
                           "--tax-mode",
                           "Creating target database Mode",
                           "Creating target database based on taxonomy of GTDB or NCBI:\n1: GTDB\t2:NCBI ",
                           typeid(int),
                           (void *) &gtdbOrNcbi,
                           "[1-2]"),
        VIRUS_TAX_ID(VIRUS_TAX_ID_ID,
                     "--virus-taxid",
                     "Taxonomy ID of virus taxon",
                     "NCBI: 10239 [Default]\nCUSTOM: Check names.dmp file ",
                     typeid(int),
                     (void *) &virusTaxId,
                     "[^[1-9]\\d*$]"),
        SEQ_MODE(SEQ_MODE_ID,
                 "--seq-mode",
                 "Taxonomy ID of virus taxon",
                 "Single-end: 1 [Default]\nPaired-end: 2\nLong read: 3",
                 typeid(int),
                 (void *) &seqMode,
                 "[1-3]"),
        MEMORY_MODE(MEMORY_MODE_ID,
                    "--memory-mode",
                    "Keeping k-mer matches in the RAM or writing into a file",
                    "Writing: 1 [Default]\n RAM:  2",
                    typeid(int),
                    (void *) &memoryMode,
                    "[1-2]"),
        REDUCED_AA(REDUCED_AA_ID,
                   "--reduce-aa",
                   "Using reduced 15 alphabets to encode amino acids. It increases sensitivity",
                   "Using 20 alphabets: 0 [Default]\n Using 15 alphabets: 1",
                   typeid(int),
                   (void *) &reducedAA,
                   "[1-2]") {
    //build_dir
    build_dir.push_back(&PARAM_THREADS);
    build_dir.push_back(&PARAM_GTDB_OR_NCBI);
    build_dir.push_back(&REDUCED_AA);

    //classify
    classify.push_back(&PARAM_THREADS);
    classify.push_back(&SEQ_MODE);
    classify.push_back(&VIRUS_TAX_ID);
    classify.push_back(&MEMORY_MODE);
    classify.push_back(&REDUCED_AA);

    //updateTargetDB

}
