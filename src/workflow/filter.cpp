#include "LocalParameters.h"
#include "QueryFilter.h"

void setFilterDefaults(LocalParameters & par){
    par.seqMode = 2;
    par.reducedAA = 0;
    par.minScore = 0.7;
    par.minCoverage = 0;
    par.minSpScore = 0;
    par.spaceMask = "11111111";
    par.hammingMargin = 0;
    par.verbosity = 3;
    par.ramUsage = 128;
    par.minCoveredPos = 4;
    par.printLog = 0;
    par.maxGap = 0;
    par.taxonomyPath = "DBDIR/taxonomy/" ;
    par.minConsCnt = 4;
    par.minConsCntEuk = 9;
    par.eukaryotaTaxId = 2759;
    par.maskMode = 0;
    par.maskProb = 0.9;
    par.matchPerKmer = 4;
    par.printMode = 1;
}

int filter(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setFilterDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    QueryFilter * queryFilter = new QueryFilter(par);
    
    queryFilter->filterReads(par);
    
    delete queryFilter;
    
    return 0;
}