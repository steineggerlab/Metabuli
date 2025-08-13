#include "IndexCreator.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "accession2taxid.h"
#include "editNames.h"
#include "fasta_validate.h"
#include "validateDatabase.h"


void setDefaults_updateDB(LocalParameters & par){
    par.makeLibrary = 0;
    par.gtdb = 0;
    par.validateInput = 0;
    par.validateDb = 0;
    par.kmerFormat = 1;
    par.skipRedundancy = 0;
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
    if (!FileUtil::directoryExists(newDbDir.c_str())) {
        FileUtil::makeDir(newDbDir.c_str());
    }

    if (par.validateInput == 1) {
        const string & fnaListFileName = par.filenames[1];
        // Read the file line by line
        ifstream file(fnaListFileName);
        if (!file.is_open()) {
            cout << "Error: Unable to open file " << fnaListFileName << endl;
            return 1;
        }
        string singleFastaName;
        while (getline(file, singleFastaName)) {
            if (!LocalUtil::isFasta(singleFastaName)) {
                cout << "Error: " << singleFastaName << " has a unsupported extension." << endl;
                cout << "       Allowed extensions are .fna, .fasta, .fa, and their gzip versions (e.g., .fna.gz)" << endl;
                file.close();
                return 1;
            }
            if (!FileUtil::fileExists(singleFastaName.c_str())) {
                cout << "Error: " << singleFastaName << " does not exist." << endl;
                file.close();
                return 1;
            }
            int validateRes = validate_fasta_file(singleFastaName.c_str(), 1);
            if (validateRes != 0) {
                cout << "Error: " << singleFastaName << " is not a valid FASTA file." << endl;
                cout << "       Please check the entries in the file" << endl;
                file.close();
                return 1;
            }
        }
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
        Debug(Debug::INFO) << "New taxonomy generated.\n";
    }

    loadDbParameters(par, oldDbDir);
    IndexCreator idxCre(par, taxonomy, par.kmerFormat);
    idxCre.setIsUpdating(true);
    idxCre.createIndex();
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
        while (fscanf(oldTaxIdListFile, "%31s", taxID) == 1) {
           taxIdSet.insert(static_cast<TaxID>(std::stoul(taxID)));
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

    idxCre.updateTaxId2SpeciesTaxId(newDbDir + "/taxID_list"); 
    idxCre.addFilesToMerge(oldDbDir + "/diffIdx", oldDbDir + "/info");
    idxCre.printFilesToMerge();
    idxCre.setMergedFileNames(newDbDir + "/diffIdx", newDbDir + "/info", newDbDir + "/split");
    idxCre.mergeTargetFiles<FilterMode::DB_CREATION>();
    
    delete taxonomy;
    cout << "Index creation completed." << endl;


    if (par.validateDb) {
        cout << "Validating the updated database..." << endl;
        if (validateDatabase(newDbDir) != 0) {
            cerr << "Database validation failed." << endl;
            return 1;
        }
        cout << "Database validation completed successfully." << endl;
    }

    return 0;
}