#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"

void setDefaults_build(LocalParameters & par){
    par.reducedAA = 0;
    par.spaceMask = "11111111";
    par.taxonomyPath = "" ;
    par.splitNum = 4096;
    par.maskProb = 0.5;
    par.maskMode = 0;
    par.bufferSize = 1'000'000'000;
}

int build(int argc, const char **argv, const Command &command){
    // Load parameters
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_build(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    string dbDirectory = par.filenames[0];
    string fastaListPath = par.filenames[1];
    string mappingFile = par.filenames[2];
    if (par.taxonomyPath.empty()) {
        par.taxonomyPath = dbDirectory + "/taxonomy/";
    } else {
        par.taxonomyPath = par.taxonomyPath + "/";
    }

    // If dbDirectory does not exist, create it
    if (!FileUtil::directoryExists(dbDirectory.c_str())) {
        FileUtil::makeDir(dbDirectory.c_str());
    }

    cout << "Taxonomy path: " << par.taxonomyPath << endl;

    IndexCreator idxCre(par, dbDirectory, fastaListPath, mappingFile);
    idxCre.createIndex(par);

    if(idxCre.getNumOfFlush() == 1) {
        cerr << "Index creation completed." << endl;
        return 0;
    }

    //Merge files
    cout << "Merge reference DB files ... " << endl;
    int numOfSplits = idxCre.getNumOfFlush();
    FileMerger merger(par);
    merger.mergeTargetFiles(par, numOfSplits);
    cerr << "Index creation completed." << endl;

    return 0;
}
