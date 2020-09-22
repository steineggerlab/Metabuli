//
// Created by KJB on 01/09/2020.
//

#ifndef ADKMER4_SEARCHER_H
#define ADKMER4_SEARCHER_H
#include "KmerExtractor.h"
#include "BitManipulateMacros.h"
#include "Mmap.h"
#include <fstream>
#include "Kmer.h"
#include "KmerExtractor.h"
#include "printBinary.h"
#include "common.h"
using namespace std;

class Searcher
{
private:
    size_t numOfSplit;
    KmerExtractor * kmerExtractor;
    size_t queryCount;
    size_t totalMatchCount;
    size_t multipleMatchCount;
    size_t perfectMatchCount;

    uint64_t lookingTarget;
    uint64_t lookingQuery;
    uint64_t nextTarget;
    uint64_t lastQuery;
    size_t lookingTargetPos;
    int isMatched;
    vector<matchedKmer> matchedKmerList;
    vector<int> closestKmers;
    int currentHamming;
    ExtractStartPoint ESP;

    uint64_t marker= ~0 & ~16777215;
    uint8_t hammingLookup[6][6]= {
            {0, 1, 1, 1, 2, 1},
            {1, 0, 1, 1, 2, 2},
            {1, 1, 0, 1, 2, 2},
            {1, 1, 1, 0, 1, 2},
            {2, 2, 2, 1, 0, 1},
            {1, 2, 2, 2, 1, 0}};
    uint8_t getHammingDistance(uint64_t kmer1, uint64_t kmer2);
    static uint64_t getNextTargetKmer(uint64_t lookingTarget, const uint16_t * targetDiffIdxList, size_t & diffIdxPos);
    void linearSearch(Kmer * kmerBuffer, size_t & bufferIdx, const MmapedData<uint16_t> & targetDiffIdxList, const MmapedData<KmerInfo> & targetInfoList, const vector<int> & taxIdList);
    void writeResultFile(vector<matchedKmer> & matchList, char * queryFileName);
public:
    void startSearch(char * queryFileName, char * targetDiffIdxFileName, char * targetInfoFileName);
    Searcher();
};
#endif //ADKMER4_SEARCHER_H
