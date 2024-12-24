#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"

void setDefaults_updateDB(LocalParameters & par){
    par.skipRedundancy = 1;
    par.reducedAA = 0;
    par.ramUsage = 128;
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

int updateDB(int argc, const char **argv, const Command &command){
    // Load parameters
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_updateDB(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    // If dbDirectory does not exist, create it
    if (!FileUtil::directoryExists(par.filenames[0].c_str())) {
        FileUtil::makeDir(par.filenames[0].c_str());
    }

    // Check if the taxonomy path exists
    if (par.taxonomyPath.empty()) {
        cerr << "Taxonomy path is not set." << endl;
        return 1;
    }
    if (!FileUtil::directoryExists(par.taxonomyPath.c_str())) {
        cerr << "Taxonomy path does not exist: " << par.taxonomyPath << endl;
        return 1;
    }
    string newDbDir = par.filenames[0];
    string oldDbDir = par.filenames[3];

    // Create index
    IndexCreator idxCre(par);
    idxCre.setIsUpdating(true);
    idxCre.createIndex(par);
    unordered_set<TaxID> taxIdSet = idxCre.getTaxIdSet();
    FILE * oldTaxIdListFile;
    if((oldTaxIdListFile = fopen((oldDbDir + "/taxID_list").c_str(),"r")) == NULL){
        cout << "Cannot open the taxID_list file of the old database" << endl;
        return 1;
    } else {
        char taxID[100];
        while(feof(oldTaxIdListFile) == 0) {
            fscanf(oldTaxIdListFile,"%s",taxID);
            taxIdSet.insert(atol(taxID));
        }
        fclose(oldTaxIdListFile);
        FILE * taxidListFile = fopen((newDbDir + "/taxID_list").c_str(), "w");
        for (auto & taxid: taxIdSet) {
            fprintf(taxidListFile, "%d\n", taxid);
        }
        fclose(taxidListFile);
    }

    // Merge index files
    cout << "Merge new and old DB files ... " << endl;
    int numOfSplits = idxCre.getNumOfFlush();
    FileMerger merger(par);
   
    for (int i = 0; i < numOfSplits; i++) {
        merger.addFilesToMerge(newDbDir + "/" + to_string(i) + "_diffIdx",
                               newDbDir + "/" + to_string(i) + "_info");
    }
    merger.addFilesToMerge(oldDbDir + "/diffIdx", oldDbDir + "/info");
    merger.setRemoveRedundancyInfo(haveRedundancyInfo(oldDbDir));
    merger.updateTaxId2SpeciesTaxId(newDbDir + "/taxID_list");
    merger.setMergedFileNames(newDbDir + "/diffIdx", newDbDir + "/info", newDbDir + "/split");
    merger.mergeTargetFiles();
    cerr << "Index creation completed." << endl;
    return 0;
}