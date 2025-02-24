#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Mmap.h"

int taxdump(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    std::string taxonomyFile = par.filenames[0];

    if (!FileUtil::fileExists(taxonomyFile.c_str())) {
        Debug(Debug::INFO) << "Taxonomy file " << taxonomyFile << " is NOT exists.\n";
        return 0;
    }

    FILE *handle = fopen(taxonomyFile.c_str(), "r");
    struct stat sb;
    if (fstat(fileno(handle), &sb) < 0) {
      Debug(Debug::ERROR) << "Failed to fstat file " << taxonomyFile << "\n";
      EXIT(EXIT_FAILURE);
    }
    char *data = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE,
                              fileno(handle), 0);
    if (data == MAP_FAILED) {
      Debug(Debug::ERROR) << "Failed to mmap file " << taxonomyFile << " with error "
                          << errno << "\n";
      EXIT(EXIT_FAILURE);
    }
    fclose(handle);
    TaxonomyWrapper *t = TaxonomyWrapper::unserialize(data);
    if (t != NULL) {
      t->setMmapData(data, sb.st_size);
    }


    t->writeNodesDmp(taxonomyFile + "_nodes.dmp");
    t->writeNamesDmp(taxonomyFile + "_names.dmp");
    t->writeMergedDmp(taxonomyFile + "_merged.dmp");

    std::cout << "Taxonomy files written." << std::endl;
    std::cout << "nodes.dmp: " << taxonomyFile << "_nodes.dmp" << std::endl;
    std::cout << "names.dmp: " << taxonomyFile << "_names.dmp" << std::endl;
    std::cout << "merged.dmp: " << taxonomyFile << "_merged.dmp" << std::endl;

    delete t;
    return 0;
}

    