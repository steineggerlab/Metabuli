#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "accession2taxid.h"
#include "editNames.h"
#include "fasta_validate.h"

void setDefaults_build(LocalParameters & par){
    par.gtdb = 0;
    par.makeLibrary = 0;
    par.reducedAA = 0;
    par.ramUsage = 128;
    par.validateInput = 0;
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
        merger.addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
                               dbDir + "/" + to_string(i) + "_info");
    }
    merger.updateTaxId2SpeciesTaxId(dbDir + "/taxID_list");
    merger.setMergedFileNames(dbDir + "/diffIdx", dbDir + "/info", dbDir + "/split");  
    merger.mergeTargetFiles();
    delete taxonomy;
    cout << "Index creation completed." << endl;
    return 0;
}
