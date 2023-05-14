#ifndef ADKMER4_SEARCHER_H
#define ADKMER4_SEARCHER_H

#include "BitManipulateMacros.h"
#include "Mmap.h"
#include <fstream>
#include "Kmer.h"
#include "SeqIterator.h"
#include "printBinary.h"
#include "common.h"
#include "NcbiTaxonomy.h"
#include "Debug.h"
#include "KmerBuffer.h"
#include "IndexCreator.h"
#include <cstdio>
#include <ctime>
#include <vector>
#include <algorithm>
#include <regex>
#include "FastSort.h"
#include "KSeqWrapper.h"
#include "LocalParameters.h"
#include <set>
#include <cmath>
#include "Match.h"
#include <unordered_set>
#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;

struct TaxonScore {
    TaxID taxId;
    float score;
    float coverage;
    int hammingDist;
    TaxonScore(TaxID taxId, float score, float coverage, int hammingDist) :
                taxId(taxId), score(score), coverage(coverage), hammingDist(hammingDist) {}
    TaxonScore() : taxId(0), score(0.0f), coverage(0.0f), hammingDist(0) {}
};

class Classifier {
protected:
    // Parameters
    int verbosity;
//    size_t localIndexBufferSize;
//    size_t localMatchBufferSize;

    string queryPath_1;
    string queryPath_2;
    string dbDir;
    string outDir;
    string jobId;

    // For spaced k-mer
    uint32_t * mask;
    uint32_t spaceNum;
    int spaceNum_int;
    int unmaskedPos[9];


    uint8_t hammingMargin;
    float minSpScore;
    int minCoveredPos;
    int maxGap;

    NcbiTaxonomy * taxonomy;
    unordered_map<TaxID, TaxID> taxId2speciesId;
    unordered_map<TaxID, TaxID> taxId2genusId;


    struct MatchBlock {
        MatchBlock(size_t start, size_t end, int id) : start(start), end(end), id(id) {}
        MatchBlock() : start(0), end(0), id(0) {}
        size_t start;
        size_t end;
        uint32_t id;
    };

    struct QueryKmerSplit {
        QueryKmerSplit(size_t start, size_t end, size_t length, const DiffIdxSplit& diffIdxSplit)
                : start(start), end(end), length(length), diffIdxSplit(diffIdxSplit) {}

        size_t start; // start idx in query k-mer list
        size_t end; // end idx in query k-mer list
        size_t length;
        DiffIdxSplit diffIdxSplit; // index in target k-mer list from where the search begins.
    };


    template<typename T>
    struct Buffer {
        T *buffer;
        size_t startIndexOfReserve;
        size_t bufferSize;

        explicit Buffer(size_t sizeOfBuffer=100) {
            buffer = (T *) malloc(sizeof(T) * sizeOfBuffer);
            bufferSize = sizeOfBuffer;
            startIndexOfReserve = 0;
        };

        size_t reserveMemory(size_t numOfKmer) {
            size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
            return offsetToWrite;
        };

        void reallocateMemory(size_t sizeOfBuffer) {
            if (sizeOfBuffer > bufferSize) {
                buffer = (T *) realloc(buffer, sizeof(T) * sizeOfBuffer);
                bufferSize = sizeOfBuffer;
            }
        };
    };

    int numOfSplit;
    unordered_map<TaxID, unsigned int> taxCounts;
    uint64_t MARKER;
    int bitsForCodon;
    uint8_t hammingLookup[8][8] = {
            {0, 1, 1, 1, 2, 1, 3, 3},
            {1, 0, 1, 1, 2, 2, 3, 2},
            {1, 1, 0, 1, 2, 2, 2, 3},
            {1, 1, 1, 0, 1, 2, 3, 3},
            {2, 2, 2, 1, 0, 1, 4, 4},
            {1, 2, 2, 2, 1, 0, 4, 4},
            {3, 3, 2, 3, 4, 4, 0, 1},
            {3, 2, 3, 3, 4, 4, 1, 0}};

    // Index reads in query file
    static void splitFASTQ(vector<SequenceBlock> & seqSegments, const string & queryPath);
    static void splitFASTA(vector<SequenceBlock> & seqSegments, const string & queryPath);

    // Extract query k-mer
    void fillQueryKmerBufferParallel(QueryKmerBuffer &kmerBuffer,
                                     MmapedData<char> &seqFile,
                                     const vector<SequenceBlock> &seqs,
                                     vector<Query> & queryList,
                                     const pair<size_t, size_t> & currentSplit,
                                     const LocalParameters &par);

    void fillQueryKmerBufferParallel(QueryKmerBuffer &kmerBuffer,
                                     const vector<SequenceBlock> &seqs,
                                     const vector<SequenceBlock> &seqs2,
                                     vector<Query> & queryList,
                                     const pair<size_t, size_t> & currentSplit,
                                     const LocalParameters &par);

    static int getMaxCoveredLength(int queryLength);

    template<typename T>
    T getQueryKmerNumber(T queryLength);

    void linearSearchParallel(
            QueryKmer *queryKmerList,
            size_t &queryKmerCnt,
            Buffer<Match> &matchBuffer,
            const LocalParameters &par);

    void compareDna(uint64_t query, vector<uint64_t> &targetKmersToCompare, vector<size_t> &selectedMatches,
                    vector<uint8_t> &selectedHammingSum, vector<uint16_t> &rightEndHammings, uint8_t frame);

    virtual uint8_t getHammingDistanceSum(uint64_t kmer1, uint64_t kmer2);

    virtual uint16_t getHammings(uint64_t kmer1, uint64_t kmer2);

    virtual uint16_t getHammings_reverse(uint64_t kmer1, uint64_t kmer2);

    void moveMatches(Match *dest, Match *src, int& matchNum);

    // Analyzing k-mer matches
    void fromMatchToClassification(const Match *matchList,
                                   size_t numOfMatches,
                                   vector<Query> & queryList,
                                   const LocalParameters &par);

    void chooseBestTaxon(uint32_t currentQuery,
                         size_t offset,
                         size_t end,
                         const Match *matchList,
                         vector<Query> & queryList,
                         const LocalParameters &par);

    void remainConsecutiveMatches(vector<const Match *> & curFrameMatches,
                                  vector<const Match *> & filteredMatches,
                                  TaxID genusId,
                                  const LocalParameters & par);

    size_t DFS(size_t curMatchIdx, const map<size_t, vector<size_t>>& linkedMatches,
             vector<size_t>& fiteredMatchIdx, size_t depth, size_t MIN_DEPTH, unordered_set<size_t>& used,
             unordered_map<size_t, size_t> & idx2depth);

    static bool isConsecutive(const Match * match1, const Match * match2);
    bool isConsecutive(const Match & match1, const Match & match2, const LocalParameters &par);



    TaxonScore getBestGenusMatches(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end,
                                   size_t offset, int queryLength, const LocalParameters &par);

    TaxonScore getBestGenusMatches(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end, size_t offset,
                                   int readLength1, int readLength2, const LocalParameters &par);

    TaxonScore getBestGenusMatches_spaced(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end, size_t offset,
                                          int readLength1, int readLength2);
    TaxonScore getBestGenusMatches_spaced(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end, size_t offset,
                                          int readLength1);

    TaxonScore scoreGenus(vector<const Match *> &filteredMatches,
                          int queryLength);

    TaxonScore scoreGenus(vector<const Match *> &filteredMatches,
                          int readLength1,
                          int readLength2);

    void scoreGenus_ExtensionScore(vector<Match> &filteredMatches,
                                   vector<vector<Match>> &matchesForEachGenus,
                                   vector<float> &scoreOfEachGenus,
                                   int readLength1, int readLength2);

    TaxonScore chooseSpecies(const std::vector<Match> &matches,
                       int queryLength,
                       vector<TaxID> &species,
                       unordered_map<TaxID, pair<int, int>> & speciesMatchRange);

    TaxonScore chooseSpecies(const std::vector<Match> &matches,
                       int read1Length,
                       int read2Length,
                       vector<TaxID> &species,
                       unordered_map<TaxID, pair<int, int>> & speciesMatchRange);

    TaxonScore scoreSpecies(const vector<Match> &matches,
                          size_t begin,
                          size_t end,
                          int queryLength);

    TaxonScore scoreSpecies(const vector<Match> &matches,
                          size_t begin,
                          size_t end,
                          int queryLength,
                          int queryLength2);

    TaxID lowerRankClassification(vector<Match> &matches, pair<int, int> &matchRange, TaxID speciesID);

    void getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> & taxCnt,
                               unordered_map<TaxID, TaxonCounts> & cladeCnt,
                               TaxID spciesID);

    TaxID BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root);

    template <typename T>
    static void loadBuffer(FILE * fp, T * buffer, size_t & bufferIdx, size_t size, int cnt){
        fseek(fp, cnt * sizeof(T), SEEK_CUR);
        fread(buffer, sizeof(T), size, fp);
        bufferIdx = 0;
    }

    template <typename T>
    static void loadBuffer(FILE * fp, T * buffer, size_t & bufferIdx, size_t size){
        fread(buffer, sizeof(T), size, fp);
        bufferIdx = 0;
    }

    // Write report
    void writeReadClassification(const vector<Query> & queryList, int queryNum, ofstream &readClassificationFile);

    void writeReportFile(const string &reportFileName, int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt);

    void writeReport(FILE *FP, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
                     unsigned long totalReads, TaxID taxID = 0, int depth = 0);

    unsigned int cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key);

    size_t AminoAcidPart(size_t kmer) const {
        return (kmer) & MARKER;
    }

    static size_t getCodonBits(size_t num) {
        return num & 0X7U;
    }

    void setMarker(uint64_t marker) {
        MARKER = marker;
        MARKER = ~MARKER;
    }

    void setNumOfBitsForCodon(int num) {
        bitsForCodon = num;
    }

    friend struct sortMatch;
public:

    void startClassify(const LocalParameters &par);

//    static uint64_t getNextTargetKmer(uint64_t lookingTarget, const uint16_t *targetDiffIdxList, size_t &diffIdxPos);

    static uint64_t getNextTargetKmer(uint64_t lookingTarget,
                                      const uint16_t * diffIdxBuffer,
                                      size_t & diffBufferIdx,
                                      size_t & totalPos);// size_t bufferSize, FILE * diffIdxFp);

    static TargetKmerInfo getKmerInfo(size_t bufferSize, FILE * kmerInfoFp, TargetKmerInfo * infoBuffer,
                              size_t & infoBufferIdx);

    explicit Classifier(LocalParameters & par);

    virtual ~Classifier();


};

struct sortMatch {
    explicit sortMatch(const Classifier * classifier) : classifier(classifier) {}
    bool operator() (const Match & a, const Match & b) const {
        if (a.qInfo.queryId < b.qInfo.queryId) return true;
        else if (a.qInfo.queryId == b.qInfo.queryId) {
            if (a.genusId < b.genusId) return true;
            else if (a.genusId == b.genusId) {
                if (a.speciesId < b.speciesId) return true;
                else if (a.speciesId == b.speciesId) {
                    if (a.qInfo.frame < b.qInfo.frame) return true;
                    else if (a.qInfo.frame == b.qInfo.frame) {
                        if (a.qInfo.position < b.qInfo.position) return true;
                        else if (a.qInfo.position == b.qInfo.position) {
                            return a.hamming < b.hamming;
//                            if (a.hamming < b.hamming) return true;
//                            else if (a.hamming == b.hamming) {
//                                if (a.rightEndHamming < b.rightEndHamming) return true;
//                                else if (a.rightEndHamming == b.rightEndHamming) {
//                                    return a.targetId < b.targetId;
//                                }
//                            }
                        }
                    }
                }
            }
        }
        return false;
    }
    const Classifier *  classifier;
};

inline uint8_t Classifier::getHammingDistanceSum(uint64_t kmer1, uint64_t kmer2) {//12345678
    uint8_t hammingSum = 0;
    hammingSum += hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)];
    hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 3U)][GET_3_BITS(kmer2 >> 3U)];
    hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 6U)][GET_3_BITS(kmer2 >> 6U)];
    hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 9U)][GET_3_BITS(kmer2 >> 9U)];
    hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 12U)][GET_3_BITS(kmer2 >> 12U)];
    hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 15U)][GET_3_BITS(kmer2 >> 15U)];
    hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 18U)][GET_3_BITS(kmer2 >> 18U)];
    hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 21U)][GET_3_BITS(kmer2 >> 21U)];
    return hammingSum;
}

inline uint16_t Classifier::getHammings(uint64_t kmer1, uint64_t kmer2) {  //hammings 87654321
    uint16_t hammings = 0;
    for (int i = 0; i < 8; i++) {
        hammings |= hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)] << 2U * i;
        kmer1 >>= bitsForCodon;
        kmer2 >>= bitsForCodon;
    }
    return hammings;
}

inline uint16_t Classifier::getHammings_reverse(uint64_t kmer1, uint64_t kmer2) {  //hammings 87654321
    uint16_t hammings = 0;
    for (int i = 0; i < 8; i++) {
        hammings |= hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)] << 2U * (7-i);
        kmer1 >>= bitsForCodon;
        kmer2 >>= bitsForCodon;
    }
    return hammings;
}

//inline uint64_t
//Classifier::getNextTargetKmer(uint64_t lookingTarget, const uint16_t *targetDiffIdxList, size_t &diffIdxPos) {
//    uint16_t fragment;
//    uint16_t check = (0x1u << 15u);
//    uint64_t diffIn64bit = 0;
//    fragment = targetDiffIdxList[diffIdxPos];
//    diffIdxPos++;
//    while (!(fragment & check)) { // 27 %
//        diffIn64bit |= fragment;
//        diffIn64bit <<= 15u;
//        fragment = targetDiffIdxList[diffIdxPos];
//        diffIdxPos++;
//    }
//    fragment &= ~check; // not; 8.47 %
//    diffIn64bit |= fragment; // or : 23.6%
//
//    return diffIn64bit + lookingTarget;
//}

inline uint64_t
Classifier::getNextTargetKmer(uint64_t lookingTarget, const uint16_t * diffIdxBuffer, size_t & diffBufferIdx, size_t & totalPos) {
//                              size_t bufferSize, FILE * diffIdxFp) {
    uint16_t fragment;
    uint16_t check = 32768; // 2^15
    uint64_t diffIn64bit = 0;
    fragment = diffIdxBuffer[diffBufferIdx++];
    totalPos ++;
    while (!(fragment & check)) { // 27 %
        diffIn64bit |= fragment;
        diffIn64bit <<= 15u;
        fragment = diffIdxBuffer[diffBufferIdx++];
        totalPos ++;
    }
    fragment &= ~check; // not; 8.47 %
    diffIn64bit |= fragment; // or : 23.6%
    return diffIn64bit + lookingTarget;
}

inline
TargetKmerInfo Classifier::getKmerInfo(size_t bufferSize, FILE * kmerInfoFp, TargetKmerInfo * infoBuffer,
                                       size_t & infoBufferIdx){
    if (unlikely(infoBufferIdx >= bufferSize)) {
        loadBuffer(kmerInfoFp, infoBuffer, infoBufferIdx, bufferSize, (int) (infoBufferIdx - bufferSize));
    }
    return infoBuffer[infoBufferIdx];
}


#endif //ADKMER4_SEARCHER_H
