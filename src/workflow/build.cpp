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
    string taxonomyDirectory = dbDirectory + "/taxonomy";
    string mappingFile = par.filenames[2];
    if (par.taxonomyPath != "") taxonomyDirectory = par.taxonomyPath;

    IndexCreator idxCre(par, dbDirectory, fastaListPath, taxonomyDirectory, mappingFile);
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