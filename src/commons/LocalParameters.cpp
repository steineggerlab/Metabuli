#include "LocalParameters.h"

LocalParameters::LocalParameters():
        Parameters(), PARAM_GTDB_OR_NCBI(PARAM_GTDB_OR_NCBI_ID, "--tax-mode", "Creating target database Mode", "Creating target database based on taxonomy of GTDB or NCBI:\n1: GTDB\t2:NCBI ",
                                         typeid(int), (void *) &gtdbOrNcbi, "[1-2]", MMseqsParameter::COMMAND_MISC)  {
    //createTargetDB
    createTargetDB.push_back(&PARAM_THREADS);
    createTargetDB.push_back(&PARAM_GTDB_OR_NCBI);

    //classify
    classify.push_back(&PARAM_THREADS);
    classify.push_back(&PARAM_GTDB_OR_NCBI);

    //updateTargetDB

}

//    :
//        Parameters()
//    {
//
//        createTargetDB.push_back(&PARAM_THREADS);
//    }