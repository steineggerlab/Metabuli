//
// Created by KJB on 03/09/2020.
//

#ifndef ADKMER4_MERGETARGETFILES_H
#define ADKMER4_MERGETARGETFILES_H
#include <vector>
#include "Mmap.h"
#include "Kmer.h"
#include <iostream>
#include "IndexCreator.h"
#include "printBinary.h"


using namespace std;

class FileMerger {
private:
    char * mergedDiffFileName;
    char * mergedInfoFileName;
    IndexCreator * cre;
public:
    FileMerger(char* mergedDiffFileName, char * mergedInfoFileNmae);
    void mergeTargetFiles(std::vector<char *> diffIdxFileNames, std::vector<char *> infoFileNames, std::vector<int> & taxIdListAtRank, std::vector<int> & taxIdList);
    void updateTargetDatabase(vector<char *> diffIdxFileNames, vector<char *> infoFileNames, vector<int> & oldTaxList, vector<int> & oldTaxAtRankList, vector<int> & newTaxListAtRank,
                              vector<int> & newTaxList, const int & seqIdOffset);
    static size_t smallest(const uint64_t *lookingKmer, const size_t &fileCnt);
    static uint64_t getNextKmer(uint64_t currentValue, const struct MmapedData<uint16_t> & diffList, size_t &idx);
};
#endif //ADKMER4_MERGETARGETFILES_H

