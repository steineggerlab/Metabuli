#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "accession2taxid.h"
#include "editNames.h"
#include "TaxonomyWrapper.h"

int createTaxDb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    std::string taxdumpDir = par.filenames[0];
    std::string outDir     = par.filenames[1];


    if (par.gtdb == 1) {
        editNames(taxdumpDir + "/names.dmp", taxdumpDir + "/taxid.map");
    }

    TaxonomyWrapper * taxonomy =  new TaxonomyWrapper(taxdumpDir + "/names.dmp",
                                                      taxdumpDir + "/nodes.dmp",
                                                      taxdumpDir + "/merged.dmp",
                                                      true);

    taxonomy->writeTaxonomyDB(outDir + "/taxonomyDB");
    delete taxonomy;
    return 0;
}