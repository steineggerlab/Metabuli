#include "LocalParameters.h"

LocalParameters::LocalParameters():
        Parameters(),
        PARAM_GTDB_OR_NCBI(PARAM_GTDB_OR_NCBI_ID,
                           "--tax-mode",
                           "Creating target database Mode",
                           "Creating target database based on taxonomy of GTDB or NCBI:\n1: GTDB\t2:NCBI ",
                           typeid(int),
                           (void *) &gtdbOrNcbi,
                           "[1-2]",
                           MMseqsParameter::COMMAND_MISC),
        VIRUS_TAX_ID(VIRUS_TAX_ID_ID,
                     "--virus-taxid",
                     "Taxonomy ID of virus taxon",
                     "NCBI: 10239 [Default]\nCUSTOM: Check names.dmp file ",
                     typeid(int),
                     (void *) &virusTaxId,
                     "[^[1-9]\\d*$]",
                     MMseqsParameter::COMMAND_COMMON),
        SEQ_MODE(SEQ_MODE_ID,
                     "--seq-mode",
                     "Taxonomy ID of virus taxon",
                     "Single-end: 1 [Default]\nPaired-end: 2\nLong read: 3",
                     typeid(int),
                     (void *) &seqMode,
                     "[1-3]",
                     MMseqsParameter::COMMAND_COMMON)
                     {
    //build_dir
    build_dir.push_back(&PARAM_THREADS);
    build_dir.push_back(&PARAM_GTDB_OR_NCBI);

    //classify
    classify.push_back(&PARAM_THREADS);
    classify.push_back(&SEQ_MODE);
    classify.push_back(&VIRUS_TAX_ID);

    //updateTargetDB

}
