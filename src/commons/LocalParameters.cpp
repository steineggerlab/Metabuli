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
                     MMseqsParameter::COMMAND_MISC) {
    //build_dir
    build_dir.push_back(&PARAM_THREADS);
    build_dir.push_back(&PARAM_GTDB_OR_NCBI);

    //classify
    classify.push_back(&PARAM_THREADS);

    //updateTargetDB

}
