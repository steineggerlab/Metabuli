#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"
#include "LocalUtil.h"
#include "fastq_info.cpp"
#include "validateDatabase.h"

void setClassifyDefaults(LocalParameters & par){
    par.syncmer = 0;
    par.smerLen = 6;
    par.kmerFormat = 1;
    // par.maxShift = 1;
    par.skipRedundancy = 0;
    par.reducedAA = 0;
    par.validateInput = 0;
    par.validateDb = 0;
    par.seqMode = 2;    
    par.minScore = 0;
    par.minSpScore = 0;
    par.hammingMargin = 0;
    par.verbosity = 3;
    par.ramUsage = 128;
    par.printLog = 0;
    par.maxGap = 0;
    par.taxonomyPath = "" ;
    par.minConsCnt = 4;
    par.minConsCntEuk = 9;
    par.maskMode = 0;
    par.maskProb = 0.9;
    par.matchPerKmer = 4;
    par.accessionLevel = 0;
    par.tieRatio = 0.95;
    par.printLineage = 0;
}

int classify(int argc, const char **argv, const Command& command) {
    LocalParameters & par = LocalParameters::getLocalInstance();
    setClassifyDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    
    string dbDir;
    if (par.seqMode == 2) {
        dbDir = par.filenames[2];
        if (FileUtil::directoryExists(par.filenames[1].c_str())) {
            cout << "Error: " << par.filenames[1] << " is a directory. Please specify a query file name." << endl;
            cout << "       For '--seq-mode 2', please provide two query files." << endl;
            exit(1);
        }

        if (!LocalUtil::isValidQueryFile(par.filenames[0])) {
            cout << "Error: " << par.filenames[0] << " is not a valid query file." << endl;
            cout << "       Allowed extensions are .fna, .fasta, .fa, .fq, .fastq, and their gzip versions (e.g., .fna.gz)" << endl;
            exit(1);
        }
        if (!LocalUtil::isValidQueryFile(par.filenames[1])) {
            cout << "Error: " << par.filenames[1] << " is not a valid query file." << endl;
            cout << "       Allowed extensions are .fna, .fasta, .fa, .fq, .fastq, and their gzip versions (e.g., .fna.gz)" << endl;
            exit(1);
        }
        if (!FileUtil::directoryExists(par.filenames[3].c_str())) {
            FileUtil::makeDir(par.filenames[3].c_str());
        }

        if (par.validateInput) {
            if (LocalUtil::isFasta(par.filenames[0])) {
                cout << "Validating FASTA file: " << par.filenames[0] << endl;
                int validateRes = validate_fasta_file(par.filenames[0].c_str(), 1);
                if (validateRes != 0) {
                    cout << "Error: " << par.filenames[0] << " is not a valid FASTA file." << endl;
                    cout << "       Please check the entries in the file" << endl;
                    return 1;
                }

                cout << "Validating FASTA file: " << par.filenames[1] << endl;
                validateRes = validate_fasta_file(par.filenames[1].c_str(), 1);
                if (validateRes != 0) {
                    cout << "Error: " << par.filenames[1] << " is not a valid FASTA file." << endl;
                    cout << "       Please check the entries in the file" << endl;
                    return 1;
                }
            }
            if (LocalUtil::isFastq(par.filenames[0])) {
                cout << "Validating FASTQ file: " << par.filenames[0] << endl;
                FASTQ_FILE* fd1=validate_single_fastq_file(par.filenames[0].c_str());
                fastq_destroy(fd1);

                cout << "Validating FASTQ file: " << par.filenames[1] << endl;
                FASTQ_FILE* fd2=validate_single_fastq_file(par.filenames[1].c_str());
                fastq_destroy(fd2);
            }
        }
        if (par.validateDb) {
            if (validateDatabase(dbDir) != 0) {
                cout << "Error: Database validation failed." << endl;
                exit(1);
            }
        }
    } else {
        dbDir = par.filenames[1];
        if (FileUtil::fileExists(par.filenames[1].c_str()) 
            && !FileUtil::directoryExists(par.filenames[1].c_str())) {
            cout << "Error: " << par.filenames[1] << " is a file. Please specify a database directory." << endl;
            cout << "       For '--seq-mode 1' and '--seq-mode 3', please provide one query file." << endl;
            exit(1);
        }
        if (!LocalUtil::isValidQueryFile(par.filenames[0])) {
            cout << "Error: " << par.filenames[0] << " is not a valid query file." << endl;
            cout << "       Allowed extensions are .fna, .fasta, .fa, .fq, .fastq, and their gzip versions (e.g., .fna.gz)" << endl;
            exit(1);
        }
        if (!FileUtil::directoryExists(par.filenames[2].c_str())) {
            FileUtil::makeDir(par.filenames[2].c_str());
        }
        if (par.validateInput) {
            if (LocalUtil::isFasta(par.filenames[0])) {
                cout << "Validating FASTA file: " << par.filenames[0] << endl;
                int validateRes = validate_fasta_file(par.filenames[0].c_str(), 1);
                if (validateRes != 0) {
                    cout << "Error: " << par.filenames[0] << " is not a valid FASTA file." << endl;
                    cout << "       Please check the entries in the file" << endl;
                    return 1;
                }
            }
            if (LocalUtil::isFastq(par.filenames[0])) {
                cout << "Validating FASTQ file: " << par.filenames[0] << endl;
                FASTQ_FILE* fd1=validate_single_fastq_file(par.filenames[0].c_str());
                fastq_destroy(fd1);
            }
        }
        if (par.validateDb) {
            if (validateDatabase(dbDir) != 0) {
                cout << "Error: Database validation failed." << endl;
                exit(1);
            }
        }
    }

    // Validate database directory
    if (!par.validateDb) {
        if (!FileUtil::directoryExists(dbDir.c_str())) {
            cout << "Error: " << dbDir << " is not found." << endl;
            exit(1);
        }
        const string & deltaIdx = dbDir + "/deltaIdx.mtbl";
        const string & diffIdx = dbDir + "/diffIdx";

        bool newFormat = false;
        if (FileUtil::fileExists(diffIdx.c_str())) {
            newFormat = false;
        } else if (FileUtil::fileExists(deltaIdx.c_str())) {
            newFormat = true;
        } else {
            cout << "Error: Neither " << diffIdx << " nor " << deltaIdx << " is found." << endl;
            cout << "       Please check the database directory." << endl;
            exit(1);
        }

        if (newFormat) {
            const string & deltaIdxSplit = dbDir + "/deltaIdxSplits.mtbl";
            if (!FileUtil::fileExists(deltaIdxSplit.c_str())) {
                cout << "Error: " << deltaIdxSplit << " is not found." << endl;
                exit(1);
            }
        } else {
            const string & info = dbDir + "/info";
            const string & split = dbDir + "/split";
            if (!FileUtil::fileExists(info.c_str())) {
                cout << "Error: " << info << " is not found." << endl;
                exit(1);
            }
            if (!FileUtil::fileExists(split.c_str())) {
                cout << "Error: " << split << " is not found." << endl;
                exit(1);
            }
        }
        const string & taxonomy = dbDir + "/taxonomyDB";
        if (!FileUtil::fileExists(taxonomy.c_str())
            && par.taxonomyPath.empty()
            && !FileUtil::fileExists((dbDir + "/taxonomy/merged.dmp").c_str())) {
            cout << "Error: taxonomy files are not found." << endl;
            cout << "       One of the followings should be provided:" << endl;
            cout << "       1. File: " << taxonomy << endl;
            cout << "       2. Dir : " << dbDir + "/taxonomy" << endl;
            cout << "       3. Specify --taxonomy-path" << endl;
            exit(1);
        }
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    Classifier * classifier = new Classifier(par);
    classifier->startClassify(par);
    delete classifier;
    return 0;
}