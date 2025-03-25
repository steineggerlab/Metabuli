#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "accession2taxid.h"
#include "editNames.h"

void setDefaults_build(LocalParameters & par){
    par.syncmer = 0;
    par.smerLen = 6;
    par.gtdb = 0;
    par.makeLibrary = 1;
    par.reducedAA = 0;
    par.ramUsage = 128;
    // par.spaceMask = "11111111";
    par.taxonomyPath = "" ;
    par.splitNum = 4096;
    par.maskProb = 0.9;
    par.maskMode = 1;
    par.accessionLevel = 0;
    // Get current date
    time_t now = time(0);
    tm *ltm = localtime(&now);
    par.dbDate = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday);
    
    // Get random alphanumeric string fore dbName from current time
    srand(time(NULL));
    string randStr = to_string(rand());
    par.dbName = randStr.substr(0, 32);
}

int build(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_build(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string & dbDir = par.filenames[0];
    if (!FileUtil::directoryExists(dbDir.c_str())) {
        FileUtil::makeDir(dbDir.c_str());
    }
    
    string taxonomyDir;
    if (par.taxonomyPath.empty()) {
        taxonomyDir = dbDir + "/taxonomy/";
    } else {
        taxonomyDir = par.taxonomyPath + "/";
    }

    if (par.gtdb == 1) {
        editNames(taxonomyDir + "/names.dmp", par.filenames[2]);
        accession2taxid(par.filenames[1], par.filenames[2]);
        par.filenames[2] = par.filenames[2].substr(0, par.filenames[2].find_last_of('.')) + ".accession2taxid";
    }

    TaxonomyWrapper * taxonomy =  new TaxonomyWrapper(taxonomyDir + "/names.dmp",
                                                      taxonomyDir + "/nodes.dmp",
                                                      taxonomyDir + "/merged.dmp",
                                                      true);

    IndexCreator idxCre(par, taxonomy);
    idxCre.createIndex(par);
    if (par.accessionLevel == 1) {
        taxonomy = idxCre.getTaxonomy();
    }
    taxonomy->writeTaxonomyDB(dbDir + "/taxonomyDB");
    
    if(idxCre.getNumOfFlush() == 1) {
        delete taxonomy;
        cout << "Index creation completed." << endl;
        return 0;
    }

    cout << "Merge reference DB files ... " << endl;
    int numOfSplits = idxCre.getNumOfFlush();
    FileMerger merger(par, taxonomy);
    for (int i = 0; i < numOfSplits; i++) {
        merger.addFilesToMerge(dbDir + "/" + to_string(i) + "_deltaIdx.mtbl", "");
    }
    merger.updateTaxId2SpeciesTaxId(dbDir + "/taxID_list");
    merger.setMergedFileNames(dbDir + "/deltaIdx.mtbl", "", dbDir + "/deltaIdxSplits.mtbl");  
    merger.mergeDeltaIdxFiles();
    delete taxonomy;
    cout << "Index creation completed." << endl;
    return 0;
}
