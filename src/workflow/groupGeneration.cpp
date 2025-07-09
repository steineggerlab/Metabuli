#include "GroupGenerator.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"

void setGroupGenerationDefaults(LocalParameters & par){
    par.reducedAA = 0;
    par.seqMode = 2;    
    par.minScore = 0;
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
    par.accessionLevel = 0;
    par.tieRatio = 0.95;

    par.thresholdK = 0.5;
    par.groupScoreThr = 0.0;
    par.voteMode = 0;
    par.majorityThr = 0.5;
    par.onlyAA = false;
}

int groupGeneration(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setGroupGenerationDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    if (par.seqMode == 2) {
        // Check if the second argument is a directory
        if (FileUtil::directoryExists(par.filenames[1].c_str())) {
            cerr << "Error: " << par.filenames[1] << " is a directory. Please specify a query file name." << endl;
            cerr << "       For '--seq-mode 2', please provide two query files." << endl;
            exit(1);
        }

        if (!FileUtil::directoryExists(par.filenames[3].c_str())) {
            FileUtil::makeDir(par.filenames[3].c_str());
        }
    } else {
        // Check if the second argument is file
        if (FileUtil::fileExists(par.filenames[1].c_str()) 
            && !FileUtil::directoryExists(par.filenames[1].c_str())) {
            cerr << "Error: " << par.filenames[1] << " is a file. Please specify a database directory." << endl;
            cerr << "       For '--seq-mode 1' and '--seq-mode 3', please provide one query file." << endl;
            exit(1);
        }

        if (!FileUtil::directoryExists(par.filenames[2].c_str())) {
            FileUtil::makeDir(par.filenames[2].c_str());
        }
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif    
    GroupGenerator * groupGenerator = new GroupGenerator(par);
    groupGenerator->startGroupGeneration(par);
    delete groupGenerator;
    return 0;
}