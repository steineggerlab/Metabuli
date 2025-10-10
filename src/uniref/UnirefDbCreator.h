#ifndef METABULI_UNIREF_DB_CREATOR_H
#define METABULI_UNIREF_DB_CREATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <ctime>

#include "LocalParameters.h"
#include "TaxonomyWrapper.h"
#include "FileUtil.h"
#include "UnirefTree.h"

// #include "tinyxml2.h"
#include "yxml.h"

struct UniRefIdx {
    int uniref90;
    int uniref50;
    UniRefIdx(int uniref90 = 0, int uniref50 = 0) : uniref90(uniref90), uniref50(uniref50) {}
    UniRefIdx() = default;
};

class UnirefDbCreator {
private:
    const LocalParameters & par;
    std::string dbDir;
    std::string unirefTaxPath;

    void parseUnirefIds(
        const std::string & unirefXmlFileName,
        std::unordered_map<std::string, std::string> & uniref100to90,
        std::unordered_map<std::string, std::string> & uniref90to50
    );

    void parseUnirefIds(
        const std::string & unirefXmlFileName,
        std::vector<std::string> &uniref100names,
        std::vector<UniRefIdx> &uniref90and50,
        std::unordered_map<std::string, int> &uniref90toIdx,
        std::unordered_map<std::string, int> &uniref50toIdx
    );

    void createUnirefDumpFiles(
        const std::unordered_map<std::string, std::string> & uniref100to90,
        const std::unordered_map<std::string, std::string> & uniref90to50
    );


    void dumpUnirefTree(
        const std::unordered_map<std::string, std::string> & uniref100to90,
        const std::unordered_map<std::string, std::string> & uniref90to50
    );

public:
    UnirefDbCreator(
        const LocalParameters &par,
        const std::string & dbDir
    );

    // void createUnirefDb();

    void createUnirefTaxonomy(const std::string & unirefXmlFileName);

    void createUnirefHierarchy(const std::string & unirefXmlFileName);

    static int getLCA(
        const std::vector<int> & ids,
        const std::unordered_map<int, std::pair<int, int>> & uniref100to90and50
    );
};

#endif
