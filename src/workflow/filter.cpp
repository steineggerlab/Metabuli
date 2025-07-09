#include "LocalParameters.h"
#include "QueryFilter.h"

void setFilterDefaults(LocalParameters & par) {
    par.reducedAA = 0;
    // par.spaceMask = "11111111";
    par.seqMode = 2;
    par.minScore = 0.5;
    par.minCoverage = 0;
    par.minSpScore = 0;
    par.hammingMargin = 0;
    par.verbosity = 3;
    par.ramUsage = 128;
    par.minCoveredPos = 4;
    par.printLog = 0;
    par.maxGap = 0;
    par.taxonomyPath = "" ;
    par.minConsCnt = 4;
    par.minConsCntEuk = 9;
    par.maskMode = 0;
    par.maskProb = 0.9;
    par.matchPerKmer = 4;
    par.printMode = 1;
    par.contamList = ""; // TODO: set default
    par.accessionLevel = 0;
    par.onlyAA = false;
}

int filter(int argc, const char **argv, const Command& command) {
    LocalParameters & par = LocalParameters::getLocalInstance();
    setFilterDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    if (par.contamList == "") {
        cerr << "Error: Contamination list is not specified." << endl;
        return 1;
    }

    QueryFilter * queryFilter = new QueryFilter(par);
    
    queryFilter->filterReads(par);
    
    delete queryFilter;
    
    return 0;
}