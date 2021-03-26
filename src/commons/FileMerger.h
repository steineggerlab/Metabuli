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
#include "common.h"


using namespace std;

class FileMerger {
private:
    char * mergedDiffFileName;
    char * mergedInfoFileName;
    char * diffIdxSplitFileName;
    IndexCreator * cre;

    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx, size_t & totalBufIdx);
    void writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx , size_t & totalBufIdx);
    void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx );

public:
    FileMerger(char* mergedDiffFileName, char * mergedInfoFileNmae, char * diffIdxSplitFileName);
    void mergeTargetFiles(std::vector<char *> diffIdxFileNames, std::vector<char *> infoFileNames, std::vector<int> & taxIdListAtRank, std::vector<int> & taxIdList);
    void updateTargetDatabase(vector<char *> diffIdxFileNames, vector<char *> infoFileNames, vector<int> & taxListAtRank, vector<int> & taxIdList, const int & seqIdOffset);
    static size_t smallest(const uint64_t *lookingKmer, const TargetKmerInfo lookingInfos[], vector<int> & taxListAtRank, const size_t &fileCnt);
    static uint64_t getNextKmer(uint64_t currentValue, const struct MmapedData<uint16_t> & diffList, size_t &idx);
};
#endif //ADKMER4_MERGETARGETFILES_H

