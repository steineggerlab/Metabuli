#ifndef METABULI_CREATE_COMMON_KMER_LIST_H
#define METABULI_CREATE_COMMON_KMER_LIST_H

#include "editNames.h"
#include "accession2taxid.h"
#include "IndexCreator.h"

void setDefaults_uniqKmers(LocalParameters & par){
    par.syncmer = 0;
    par.smerLen = 5;
    par.reducedAA = 0;
    par.ramUsage = 128;
    par.validateInput = 0;
    par.taxonomyPath = "" ;
    par.splitNum = 4096;
    time_t now = time(0);
    tm *ltm = localtime(&now);
    par.dbDate = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday);
    srand(time(NULL));
    string randStr = to_string(rand());
    par.dbName = randStr.substr(0, 32);
}

int create_unique_kmer_list(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_uniqKmers(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string & dbDir = par.filenames[0];

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

    cout << "Creating unique k-mer list in " << dbDir << endl;
    IndexCreator idxCre(par, 4);
    idxCre.createUniqueKmerIndex();  
    return 0;
}

#endif // METABULI_CREATE_COMMON_KMER_LIST_H