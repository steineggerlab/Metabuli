#include "LocalParameters.h"
#include "UnirefDbCreator.h"
#include "IndexCreator.h"
#include "validateDatabase.h"

void setDefaults_create_unirefdb(LocalParameters & par){
    par.ramUsage = 128;
    par.validateDb = 0;
    par.taxonomyPath = "";
    par.splitNum = 4096;
    par.unirefXml = "";
}

int create_unirefdb(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_create_unirefdb(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string & dbDir        = par.filenames[0];
    const std::string & unirefFasta  = par.filenames[1];
    const std::string & unirefTreeDb = par.filenames[2];

    if (!FileUtil::fileExists(unirefFasta.c_str())) {
        std::cerr << "Error: UniRef100 FASTA file " << unirefFasta << " does not exist." << std::endl;
        return 1;
    }
    if (!FileUtil::fileExists(unirefTreeDb.c_str())) {
        std::cerr << "Error: UniRef tree database " << unirefTreeDb << " does not exist." << std::endl;
        return 1;
    }
    if (!FileUtil::directoryExists(dbDir.c_str())) {
        if (!FileUtil::makeDir(dbDir.c_str())) {
            std::cerr << "Error: Could not create output directory " << dbDir << std::endl;
            return 1;
        }
    }

    time_t start = time(nullptr);
    std::cout << "Loading Taxonomy" << std::endl;
    UnirefTree * unirefTree = UnirefTree::openUnirefTree(unirefTreeDb);
    std::cout << "Taxonomy loaded in " << time(nullptr) - start << " seconds." << std::endl;

    IndexCreator idxCre(par, unirefTree, 4);
    idxCre.createLcaKmerIndex();

    if (par.validateDb)
    {
        std::cout << "Validating the created database..." << std::endl;
        if (validateDatabase(dbDir) != 0)
        {
            std::cerr << "Database validation failed." << std::endl;
            return 1;
        }
        std::cout << "Database validation completed successfully." << std::endl;
    }

    return 0;
}
