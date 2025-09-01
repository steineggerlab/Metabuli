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

    const std::string & dbDir = par.filenames[0];
    std::string taxonomyDir;
    
    if (par.taxonomyPath.empty()) {
        if (par.unirefXml.empty()) {
            std::cerr << "Error: --uniref-xml option is required when taxonomy path is not provided." << std::endl;
            return 1;
        }
        UnirefDbCreator unirefDbCreator(par, dbDir);
        unirefDbCreator.createUnirefTaxonomy(par.unirefXml);
        taxonomyDir = dbDir + "/taxonomy";
    } else {
        taxonomyDir = par.taxonomyPath;
    }

    return 0;

    time_t start = time(nullptr);
    std::cout << "Loading Taxonomy" << std::endl;
    TaxonomyWrapper * taxonomy =  new TaxonomyWrapper(taxonomyDir + "/names.dmp",
                                                      taxonomyDir + "/nodes.dmp",
                                                      taxonomyDir + "/merged.dmp",
                                                      false);
    std::cout << "Taxonomy loaded in " << time(nullptr) - start << " seconds." << std::endl;

    IndexCreator idxCre(par, taxonomy, 4);
    idxCre.createLcaKmerIndex();
    taxonomy->writeTaxonomyDB(dbDir + "/taxonomyDB");
    delete taxonomy;

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
