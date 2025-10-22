#include "UnirefClassifier.h"
#include "fasta_validate.h"

void setDefaults_uniref(LocalParameters & par){
    par.ramUsage = 128;
    par.matchPerKmer = 3;
    par.validateInput = 0;
}

int assign_uniref(int argc, const char **argv, const Command& command) {
    LocalParameters & par = LocalParameters::getLocalInstance();
    setDefaults_uniref(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    const string & queryFile  = par.filenames[0];
    const string & dbDir      = par.filenames[1];
    const string & unirefTree = par.filenames[2];
    const string & outDir     = par.filenames[3];

    if (!FileUtil::directoryExists(outDir.c_str())) {
        FileUtil::makeDir(outDir.c_str());
    }

    if (par.validateInput) {
        if (LocalUtil::isFasta(queryFile)) {
            cout << "Validating FASTA file: " << queryFile << endl;
            int validateRes = validate_fasta_file(queryFile.c_str(), 1);
            if (validateRes != 0) {
                cout << "Error: " << queryFile << " is not a valid FASTA file." << endl;
                cout << "       Please check the entries in the file" << endl;
                return 1;
            }
        } else {
            std::cout << "Error: Currently only FASTA format is supported for UniRef classification." << std::endl;
            return 1;
        }
    }

    UnirefTree * tree = UnirefTree::openUnirefTree(unirefTree);
    UnirefClassifier * classifier = new UnirefClassifier(par, tree);
    classifier->classify();
    delete classifier;
    return 0;
}