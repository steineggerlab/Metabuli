#include "LocalParameters.h"
#include "FileUtil.h"
#include "TaxonomyWrapper.h"

#include <iostream>
#include <unordered_map>
#include <unordered_set>

int uniref2taxonomy(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    string mappingFileName = par.filenames[0];
    string dbDir = par.filenames[1];

    if (!FileUtil::fileExists(mappingFileName)) {
        std::cerr << "Error: Mapping file " << mappingFileName << " does not exist." << std::endl;
        return 1;
    }

    if (!FileUtil::directoryExists(dbDir)) {
        std::cerr << "Error: Database directory " << dbDir << " does not exist." << std::endl;
        return 1;
    }

    // Open the mapping file
    std::ifstream mappingFile(mappingFileName);
    if (!mappingFile.is_open()) {
        std::cerr << "Error: Could not open mapping file " << mappingFileName << std::endl;
        return 1;
    }

    std::unordered
    std::unordered_map<std::string, std::string> uniref100to90;
    std::unordered_map<std::string, std::string> uniref90to50;
    // Read the mapping file and process it
    std::string line;
    size_t id = 1;

    while (std::getline(mappingFile, line)) {
        std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
        const std::string &uniref100 = columns[7];
        const std::string &uniref90 = columns[8];
        const std::string &uniref50 = columns[9];

        uniref100to90[uniref100] = uniref90;
        uniref90to50[uniref90] = uniref50;

   }
    return 0;
}