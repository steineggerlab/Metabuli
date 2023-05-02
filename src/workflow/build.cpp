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
    if (par.tinfoPath.empty()) {
        par.tinfoPath = dbDirectory + "/prodigal/";
    } else {
        par.tinfoPath = par.tinfoPath + "/";
    }
    cout << par.maskProb << endl;
    
    // If the prodigal directory does not exist, create it
    if (!FileUtil::directoryExists(par.tinfoPath.c_str())) {
        FileUtil::makeDir(par.tinfoPath.c_str());
    }
    cout << "Taxonomy path: " << par.taxonomyPath << endl;
    cout << "Tinfo path: " << par.tinfoPath << endl;

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
