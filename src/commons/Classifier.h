//
// Created by KJB on 01/09/2020.
//

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

#include <time.h>
#include <vector>
#include <algorithm>
#include <regex>
#include "FastSort.h"
#include "KSeqWrapper.h"
#include "LocalParameters.h"
#include <set>
#include <cmath>
#include "Match.h"


#define AminoAcid(x) (size_t)((x) & (~0 & ~16777215))
using namespace std;

struct Counts {
    int classificationCnt;
    int correct;
    int highRank;

    //number of targets at each rank
    int subspeciesTargetNumber;
    int speciesTargetNumber;
    int genusTargetNumber;
    int familyTargetNumber;
    int orderTargetNumber;
    int classTargetNumber;
    int phylumTargetNumber;
    int superkingdomTargetNumber;

    //number of classification at each rank
    int subspeciesCnt_try;
    int speciesCnt_try;
    int genusCnt_try;
    int familyCnt_try;
    int orderCnt_try;
    int classCnt_try;
    int phylumCnt_try;
    int superkingdomCnt_try;


    //number of correct classifications at each rank
    int subspeciesCnt_correct;
    int speciesCnt_correct;
    int genusCnt_correct;
    int familyCnt_correct;
    int orderCnt_correct;
    int classCnt_correct;
    int phylumCnt_correct;
    int superkingdomCnt_correct;
};


class Classifier {
private:

    struct ScrCov {
        float score;
        float coverage;

        ScrCov(float score, float coverage) : score(score), coverage(coverage) {}

        ScrCov() : score(0.f), coverage(0.f) {}
    };

    struct QueryInfo {
        int queryId;
        bool isClassified;
        string name;
        int taxId;
        float coverage;
        unordered_map<TaxID, int> taxCnt; ///how about using it for REPORTFILE? --> k-mer count
        size_t queryLength;

        QueryInfo(int queryId, bool isClassified, string name, int taxId, float coverage, size_t queryLength)
                : queryId(queryId), isClassified(isClassified), name(name), taxId(taxId), coverage(coverage),
                  queryLength(queryLength) {}

        QueryInfo() {}

        bool operator==(const int Id) const {
            if (Id == queryId)
                return true;
            return false;
        }
    };


    struct MatchBlock {
        MatchBlock(size_t start, size_t end, int id) : start(start), end(end), id(id) {}

        MatchBlock() : start(0), end(0), id(0) {}

        size_t start;
        size_t end;
        uint32_t id;
    };

    struct QueryKmerSplit {
        QueryKmerSplit(size_t start, size_t end, size_t length, DiffIdxSplit diffIdxSplit)
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

        explicit Buffer(size_t sizeOfBuffer) {
            buffer = (T *) malloc(sizeof(T) * sizeOfBuffer);
            bufferSize = sizeOfBuffer;
            startIndexOfReserve = 0;
        };

        size_t reserveMemory(size_t numOfKmer) {
            size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
            return offsetToWrite;
        };
    };

    int numOfSplit;
    //SeqIterator * seqIterator;
    size_t queryCount;
    size_t perfectMatchCount;
    size_t selectedMatchCount;

    // performance test
    Counts counts;

    size_t subspCnt;
    size_t speciesCnt;
    size_t genusCnt;
    size_t familyCnt;
    size_t orderCnt;
    size_t classCnt;
    size_t phylumCnt;
    size_t superCnt;

    size_t correctCnt;
    size_t perfectCnt;
    size_t classifiedCnt;

    vector <size_t> closestKmers;

    vector <QueryInfo> queryInfos;
    unordered_map<TaxID, unsigned int> taxCounts;

    // Variables changed by the number of alphabets
    int reducedAA;
    uint64_t MARKER;
    int bitsForCodon;

    uint8_t hammingLookup[11][11] = {
            {0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3},
            {1, 0, 1, 1, 2, 1, 2, 2, 2, 2, 3},
            {1, 1, 0, 1, 2, 2, 1, 2, 2, 3, 2},
            {1, 1, 1, 0, 2, 2, 2, 1, 1, 3, 3},
            {1, 2, 2, 2, 0, 1, 1, 1, 2, 4, 4},
            {2, 1, 2, 2, 1, 0, 1, 2, 4, 4, 4},
            {2, 2, 1, 2, 1, 1, 0, 2, 4, 4, 4},
            {2, 2, 2, 1, 1, 2, 2, 0, 1, 4, 4},
            {2, 2, 2, 1, 2, 4, 4, 1, 0, 4, 4},
            {3, 2, 3, 3, 4, 4, 4, 4, 4, 0, 4},
            {3, 3, 2, 3, 4, 4, 4, 4, 4, 4, 0}};

    // Extract query k-mer
    void fillQueryKmerBufferParallel(QueryKmerBuffer &kmerBuffer, MmapedData<char> &seqFile, vector <Sequence> &seqs,
                                     bool *checker, size_t &processedSeqCnt, Query *queryList,
                                     const LocalParameters &par);

    void fillQueryKmerBufferParallel_paired(QueryKmerBuffer &kmerBuffer,
                                            MmapedData<char> &seqFile1,
                                            MmapedData<char> &seqFile2,
                                            vector <Sequence> &seqs,
                                            vector <Sequence> &seqs2,
                                            bool *checker,
                                            size_t &processedSeqCnt,
                                            Query *queryList,
                                            size_t numOfSeq,
                                            const LocalParameters &par);

    static int getMaxCoveredLength(int queryLength);

    static int getQueryKmerNumber(int queryLength);


    // Linear search
    static bool compareForLinearSearch(const QueryKmer &a, const QueryKmer &b);

    void linearSearchParallel(
            QueryKmer *queryKmerList,
            size_t &queryKmerCnt,
            const char *targetDiffIdxList,
            const char *targetInfoList,
            const char *diffIdxSplits,
            Buffer<Match> &matchBuffer,
            const vector<int> &taxIdList,
            const vector<int> &speciesTaxIdList,
            const vector <TaxID> &genusTaxIdList,
            FILE *matchFile,
            const LocalParameters &par);

    void linearSearchParallel2(
            QueryKmer *queryKmerList,
            size_t &queryKmerCnt,
            const char *targetDiffIdxList,
            const char *targetInfoList,
            const char *diffIdxSplits,
            Buffer<Match> &matchBuffer,
            const vector<int> &taxIdList,
            const vector<int> &speciesTaxIdList,
            const vector <TaxID> &genusTaxIdList,
            FILE *matchFile,
            const LocalParameters &par);

    void compareDna(uint64_t query,
                    vector <uint64_t> &targetKmersToCompare, size_t startIdx,
                    vector <size_t> &selectedMatches, vector <uint8_t> &selectedHammingSum,
                    vector <uint16_t> &rightEndHammings, int i);

    uint8_t getHammingDistanceSum(uint64_t kmer1, uint64_t kmer2);

    uint16_t getHammings(uint64_t kmer1, uint64_t kmer2);

    void moveMatches(Match *dest, Match *src, int &matchNum);

    // Analyzing k-mer matches
    void analyseResultParallel(NcbiTaxonomy &ncbiTaxonomy,
                               Match *matchList,
                               size_t numOfMatches,
                               int seqNum,
                               Query *queryList,
                               const LocalParameters &par);

    static bool sortByGenusAndSpecies2(const Match &a, const Match &b);

    static bool sortBySpecies(const Match &a, const Match &b);

    void chooseBestTaxon(NcbiTaxonomy &ncbiTaxonomy,
                         uint32_t currentQuery,
                         size_t offset,
                         size_t end,
                         Match *matchList,
                         Query *queryList,
                         const LocalParameters &par);

    void chooseBestTaxon_index(NcbiTaxonomy &ncbiTaxonomy,
                               uint32_t currentQuery,
                               size_t offset,
                               size_t end,
                               Match *matchList,
                               Query *queryList,
                               const LocalParameters &par);

    static int getMatchesOfTheBestGenus(vector <Match> &matchesForMajorityLCA, Match *matchList, size_t end,
                                        size_t offset, int queryLength, float &bestScore);

    static int getMatchesOfTheBestGenus_index(vector <size_t> &matchesForMajorityLCA, Match *matchList, size_t end,
                                              size_t offset, int queryLength, float &bestScore);

    static int getMatchesOfTheBestGenus_paired(vector <Match> &matchesForMajorityLCA, Match *matchList, size_t end,
                                               size_t offset, int readLength1, int readLength2, float &bestScore);

    static void constructMatchCombination(vector <Match> &filteredMatches,
                                          vector <vector<Match>> &matchesForEachGenus,
                                          vector<float> &scoreOfEachGenus,
                                          int queryLength);

    static void constructMatchCombination_index(vector <size_t> &filteredMatchesIndex,
                                                Match *matchList,
                                                vector <vector<size_t>> &matchesForEachGenus,
                                                vector<float> &scoreOfEachGenus,
                                                int queryLength);

    static void constructMatchCombination_paired(vector <Match> &filteredMatches,
                                                 vector <vector<Match>> &matchesForEachGenus,
                                                 vector<float> &scoreOfEachGenus,
                                                 int readLength1, int readLength2);

    static bool sortMatchesByPos(const Match &a, const Match &b);


    static TaxID classifyFurther2(const std::vector <Match> &matches,
                                  NcbiTaxonomy &taxonomy,
                                  float maxKmerCnt);

    static void chooseSpecies(const std::vector <Match> &matches,
                              NcbiTaxonomy &taxonomy,
                              int queryLength,
                              ScrCov &speciesScrCov,
                              vector <TaxID> &species);

    static void classifyFurther_paired(const std::vector <Match> &matches,
                                       NcbiTaxonomy &taxonomy,
                                       int read1Length,
                                       int read2Length,
                                       ScrCov &speciesScrCov,
                                       vector <TaxID> &species);

    static ScrCov scoreTaxon(const vector <Match> &matches,
                             size_t begin,
                             size_t end,
                             int queryLength);

    static ScrCov scoreTaxon_paired(const vector <Match> &matches,
                                    size_t begin,
                                    size_t end,
                                    int queryLength,
                                    int queryLength2);

    // Write report
    void writeReadClassification(Query *queryList, int queryNum, ofstream &readClassificationFile);

    void writeReportFile(const string &reportFileName, NcbiTaxonomy &ncbiTaxonomy, int numOfQuery);

    void writeReport(FILE *fp, const NcbiTaxonomy &ncbiTaxonomy, const unordered_map <TaxID, TaxonCounts> &cladeCounts,
                     unsigned long totalReads, TaxID taxID = 0, int depth = 0);

    void writeMatches(Buffer<Match> &matchBuffer, FILE *matchFile);

    static bool compareForWritingMatches(const Match &a, const Match &b);

    unsigned int cladeCountVal(const std::unordered_map <TaxID, TaxonCounts> &map, TaxID key);

    size_t aminoAcid(size_t kmer) {
        return (kmer) & MARKER;
    }

    size_t getThreeBits(size_t num){
        return num & 0X7U;
    }
public:
    void startClassify(const char *queryFileName, const char *targetDiffIdxFileName, const char *targetInfoFileName,
                       const char *diffIdxSplitFileName, vector<int> &taxIdList, const LocalParameters &par,
                       NcbiTaxonomy &taxonomy);

    static uint64_t getNextTargetKmer(uint64_t lookingTarget, const uint16_t *targetDiffIdxList, size_t &diffIdxPos);

    int getNumOfSplits() const;

    Classifier(LocalParameters & par);

    ~Classifier();

    void performanceTest(NcbiTaxonomy &ncbiTaxonomy, Query *queryList, int numOfquery, vector<int> &wrongs);

    void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy &ncbiTaxonomy, vector<int> &wrongs, int i);
};

inline uint8_t Classifier::getHammingDistanceSum(uint64_t kmer1, uint64_t kmer2) {//12345678
    uint8_t hammingSum = 0;
    if(reducedAA == 0) {
        hammingSum += hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)];
        hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 3U)][GET_3_BITS(kmer2 >> 3U)];
        hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 6U)][GET_3_BITS(kmer2 >> 6U)];
        hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 9U)][GET_3_BITS(kmer2 >> 9U)];
        hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 12U)][GET_3_BITS(kmer2 >> 12U)];
        hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 15U)][GET_3_BITS(kmer2 >> 15U)];
        hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 18U)][GET_3_BITS(kmer2 >> 18U)];
        hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 21U)][GET_3_BITS(kmer2 >> 21U)];
    } else {
        hammingSum += hammingLookup[GET_4_BITS(kmer1)][GET_4_BITS(kmer2)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 4U)][GET_4_BITS(kmer2 >> 4U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 8U)][GET_4_BITS(kmer2 >> 8U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 12U)][GET_4_BITS(kmer2 >> 12U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 16U)][GET_4_BITS(kmer2 >> 16U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 20U)][GET_4_BITS(kmer2 >> 20U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 24U)][GET_4_BITS(kmer2 >> 24U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 28U)][GET_4_BITS(kmer2 >> 28U)];
    }
    return hammingSum;
}

inline uint16_t Classifier::getHammings(uint64_t kmer1, uint64_t kmer2) {  //hammings 87654321
    uint16_t hammings = 0;
    if(reducedAA == 1){
        for (int i = 0; i < 8; i++) {
            hammings |= hammingLookup[GET_4_BITS(kmer1)][GET_4_BITS(kmer2)] << 2U * i;
            kmer1 >>= 4U;
            kmer2 >>= 4U;
        }
    }
    for (int i = 0; i < 8; i++) {
        hammings |= hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)] << 2U * i;
        kmer1 >>= bitsForCodon;
        kmer2 >>= bitsForCodon;
    }
    return hammings;
}

inline uint64_t
Classifier::getNextTargetKmer(uint64_t lookingTarget, const uint16_t *targetDiffIdxList, size_t &diffIdxPos) {
    uint16_t fragment;
    uint16_t check = (0x1u << 15u);
    uint64_t diffIn64bit = 0;

    fragment = targetDiffIdxList[diffIdxPos];
    diffIdxPos++;
    while (!(fragment & check)) { // 27 %
        diffIn64bit |= fragment;
        diffIn64bit <<= 15u;
        fragment = targetDiffIdxList[diffIdxPos];
        diffIdxPos++;
    }
    fragment &= ~check; // not; 8.47 %
    diffIn64bit |= fragment; // or : 23.6%

    return diffIn64bit + lookingTarget;
}

#endif //ADKMER4_SEARCHER_H
