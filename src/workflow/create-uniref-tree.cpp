#include "LocalParameters.h"
#include "UnirefTree.h"
#include "Util.h"
#include "FileUtil.h"

void setDefaults_create_uniref_tree(LocalParameters & par){
}

int create_uniref_tree(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_create_uniref_tree(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string & dbDir = par.filenames[0];
    const std::string & xmlFileName = par.filenames[1];

    if (!FileUtil::fileExists(xmlFileName.c_str())) {
        std::cerr << "Error: UniRef XML file " << xmlFileName << " does not exist." << std::endl;
        return 1;
    }

    if (!FileUtil::directoryExists(dbDir.c_str())) {
        if (!FileUtil::makeDir(dbDir.c_str())) {
            std::cerr << "Error: Could not create output directory " << dbDir << std::endl;
            return 1;
        }
    }

    if (par.unirefNumbers.empty()) {
        std::cerr << "Error: --uniref-size option is required." << std::endl;
        return 1;
    }

    std::vector<std::string> sizes = Util::split(par.unirefNumbers, ",");
    if (sizes.size() != 3 || !Util::isNumber(sizes[0]) || !Util::isNumber(sizes[1]) || !Util::isNumber(sizes[2])) {
        std::cerr << "Error: --uniref-size must be in the format of <num100>,<num90>,<num50>." << std::endl;
        return 1;
    }
    uint32_t uniref100num = std::stoul(sizes[0]);
    uint32_t uniref90num  = std::stoul(sizes[1]);
    uint32_t uniref50num  = std::stoul(sizes[2]);

    UnirefTree tree(xmlFileName, uniref100num, uniref90num, uniref50num);
    // std::string dumpName = dbDir + "/uniref_tree.dmp";
    // tree.dumpUnirefTree(dumpName);
    std::string idxFileName = dbDir + "/uniref_tree.mtbl";
    tree.writeUnirefTree(idxFileName);
    return 0;
}