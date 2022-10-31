#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>

void setDefaults_build(LocalParameters & par){
    par.reducedAA = 0;
    par.spaceMask = "11111111";
    par.taxonomyPath = "" ;
}

int build(int argc, const char **argv, const Command &command){
    // Load parameters
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_build(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    string dbDirectory = par.filenames[0];
    string fastaListPath = par.filenames[1];
    string mappingFile = par.filenames[2];
    if (par.taxonomyPath.empty()) par.taxonomyPath = dbDirectory + "/taxonomy/";
    if (par.tinfoPath.empty()) par.tinfoPath = dbDirectory + "/prodigal/";
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
//    int numOfSplits = idxCre.getNumOfFlush();
    FileMerger merger(par);
    merger.mergeTargetFiles(par, 67);
    cerr << "Index creation completed." << endl;
    return 0;
}