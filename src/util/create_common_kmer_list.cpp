#ifndef METABULI_CREATE_COMMON_KMER_LIST_H
#define METABULI_CREATE_COMMON_KMER_LIST_H

#include "editNames.h"
#include "accession2taxid.h"
#include "IndexCreator.h"


void setDefaults(LocalParameters & par){
    par.kmerFormat = 3;
    par.syncmer = 0;
    par.smerLen = 5;
    par.gtdb = 0;
    par.makeLibrary = 0;
    par.reducedAA = 0;
    par.ramUsage = 128;
    par.validateInput = 0;
    par.validateDb = 0;
    par.taxonomyPath = "" ;
    par.splitNum = 4096;
    par.maskProb = 0.9;
    par.maskMode = 0;
    par.accessionLevel = 0;
    time_t now = time(0);
    tm *ltm = localtime(&now);
    par.dbDate = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday);
    srand(time(NULL));
    string randStr = to_string(rand());
    par.dbName = randStr.substr(0, 32);
}

int create_common_kmer_list(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string & dbDir = par.filenames[0];
    const string taxonomyDir = par.filenames[3];

    #ifdef OPENMP
        omp_set_num_threads(par.threads);
    #endif

    if (!FileUtil::directoryExists(dbDir.c_str())) {
        FileUtil::makeDir(dbDir.c_str());
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

    TaxonomyWrapper * taxonomy =  
        new TaxonomyWrapper(taxonomyDir + "/names.dmp",
                            taxonomyDir + "/nodes.dmp",
                            taxonomyDir + "/merged.dmp",
                            true);

    IndexCreator idxCre(par, taxonomy, par.kmerFormat);
    idxCre.createCommonKmerIndex();  
    delete taxonomy;
    return 0;
}

#endif // METABULI_CREATE_COMMON_KMER_LIST_H