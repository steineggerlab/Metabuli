#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "accession2taxid.h"
#include "editNames.h"


void setDefaults_updateDB(LocalParameters & par){
    par.makeLibrary = 1;
    // par.skipRedundancy = 1;
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
    string newDbDir = par.filenames[0];
    string oldDbDir = par.filenames[3];

    // If dbDirectory does not exist, create it
    if (!FileUtil::directoryExists(newDbDir.c_str())) {
        FileUtil::makeDir(newDbDir.c_str());
    }

    if (par.gtdb == 1) {
        accession2taxid(par.filenames[1], par.filenames[2]);
        par.filenames[2] = par.filenames[2].substr(0, par.filenames[2].find_last_of('.')) + ".accession2taxid";
    }
    
    // Load older taxonomy DB
    Debug(Debug::INFO) << "Loading taxonomy DB from " << oldDbDir << " ... ";
    TaxonomyWrapper * taxonomy = loadTaxonomy(oldDbDir);
    Debug(Debug::INFO) << "done.\n";
    FileUtil::copyFile(oldDbDir + "/acc2taxid.map", newDbDir + "/acc2taxid.map");

    // Make a new taxonomy DB if new taxa are added
    if (!par.newTaxa.empty()) {
        Debug(Debug::INFO) << "Adding new taxa to the taxonomy DB.\n";
        taxonomy->checkNewTaxa(par.newTaxa);
        std::vector<NewTaxon> newTaxaList;
        TaxonomyWrapper::getListOfTaxa(par.newTaxa, newTaxaList);
        TaxonomyWrapper * newTaxonomy = taxonomy->addNewTaxa(newTaxaList);
        delete taxonomy;
        taxonomy = newTaxonomy;
        // taxonomy->writeNamesDmp(newDbDir + "/newnodes.dmp");
        Debug(Debug::INFO) << "New taxonomy generated.\n";
    }

    IndexCreator idxCre(par, taxonomy);
    idxCre.setIsUpdating(true);
    idxCre.createIndex(par);
    if (par.accessionLevel == 1) {
        taxonomy = idxCre.getTaxonomy();
    }

    if (taxonomy->IsExternalData()) {
        FileUtil::copyFile(oldDbDir + "/taxonomyDB", newDbDir + "/taxonomyDB");
    } else {
        taxonomy->writeTaxonomyDB(newDbDir + "/taxonomyDB");
    }
    
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
    cout << "Merge new and old DB files " << endl;;
    int numOfSplits = idxCre.getNumOfFlush();
    FileMerger merger(par, taxonomy);
    for (int i = 0; i < numOfSplits; i++) {
        merger.addFilesToMerge(newDbDir + "/" + to_string(i) + "_diffIdx",
                               newDbDir + "/" + to_string(i) + "_info");
    }
    merger.addFilesToMerge(oldDbDir + "/diffIdx", oldDbDir + "/info");
    merger.setMergedFileNames(newDbDir + "/diffIdx", newDbDir + "/info", newDbDir + "/split");
    merger.printFilesToMerge();
    merger.setRemoveRedundancyInfo(haveRedundancyInfo(oldDbDir));
    merger.updateTaxId2SpeciesTaxId(newDbDir + "/taxID_list");
    Debug(Debug::INFO) << "Species-level taxonomy IDs are prepared.\n";
    
    merger.mergeTargetFiles();
    delete taxonomy;
    cout << "Index creation completed." << endl;
    return 0;
}