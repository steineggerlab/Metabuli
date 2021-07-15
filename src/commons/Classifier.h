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


#define AminoAcid(x) (size_t)((x) & (~0 & ~16777215))
using namespace std;

class Classifier
{
private:
    struct taxNode {
        void set(const double weightInput, const bool isCandidateInput, const TaxID & childTaxonInput) {
            weight = weightInput;
            isCandidate = isCandidateInput;
            childTaxon = childTaxonInput;
        }

        void update(const double weightToAdd, const TaxID & childTaxonInput) {
            if (childTaxon != childTaxonInput) { //isCandidate가 뭐야??
                isCandidate = true;
                childTaxon = childTaxonInput;
            }
            weight += weightToAdd;
        }

        // these will be filled when iterating over all contributing lineages
        double weight;
        bool isCandidate;
        TaxID childTaxon;
    };

    typedef struct ConsecutiveMatches{
        ConsecutiveMatches(uint32_t begin, uint32_t end, int matchCnt_, int hamming,
                           uint32_t gapCnt, size_t bi, size_t ei, uint8_t frame_, int score_ = 0)
            : begin(begin), end(end), matchCnt(matchCnt_), hamming(hamming), diffPosCnt(gapCnt), beginIdx(bi), endIdx(ei), frame(frame_), score(score_) {}
        uint32_t begin; //start position on query sequence
        uint32_t end; //end position on query sequence
        int matchCnt;
        int hamming; //hamming sum
        int diffPosCnt; //gap sum
        size_t beginIdx; //beginning index on matchList
        size_t endIdx; //end index
        int score;
        uint8_t frame;
    }ConsecutiveMatches;

    struct QueryInfo{
        int queryId;
        bool isClassified;
        string name;
        int taxId;
        float coverage;
        unordered_map<TaxID,int> taxCnt; ///how about using it for REPORTFILE? --> k-mer count
        size_t queryLength;
        QueryInfo(int queryId, bool isClassified, string name, int taxId, float coverage, size_t queryLength)
        : queryId(queryId), isClassified(isClassified), name(name), taxId(taxId), coverage(coverage), queryLength(queryLength) {}
        QueryInfo(){}
        bool operator == (const int Id) const{
            if(Id == queryId)
                return true;
            return false;
        }
    };

    struct Match{ //16byte
        uint32_t queryId;
        int taxID;
        int genusTaxID;
        uint16_t position;
        uint8_t frame;
        uint8_t hamming;
        int red; ///TODO remove it later
    };

    struct QueryKmerSplit{
        QueryKmerSplit(size_t start, size_t end, size_t length, DiffIdxSplit diffIdxSplit)
            : start(start), end(end), length(length), diffIdxSplit(diffIdxSplit) { }
        size_t start; // start idx in query k-mer list
        size_t end; // end idx in query k-mer list
        size_t length;
        DiffIdxSplit diffIdxSplit; // index in target k-mer list from where the search begins.
    };

    template <typename T>
    struct Buffer{
        T * buffer;
        size_t startIndexOfReserve;
        size_t bufferSize;

        Buffer(size_t sizeOfBuffer){
            buffer = (T *) malloc(sizeof(T) * sizeOfBuffer);
            bufferSize = sizeOfBuffer;
            startIndexOfReserve = 0;
        };

        size_t reserveMemory(size_t numOfKmer){
            size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
            return offsetToWrite;
        };
    };

    struct MatchBuffer{
        Match * buffer;
        size_t startIndexOfReserve;
        size_t bufferSize;

        MatchBuffer(size_t sizeOfBuffer){
            buffer = (Match *) malloc(sizeof(Match) * sizeOfBuffer);
            bufferSize = sizeOfBuffer;
            startIndexOfReserve = 0;
        };

        size_t reserveMemory(size_t numOfKmer){
            size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
            return offsetToWrite;
        };
    };


    int numOfSplit;
    SeqIterator * seqIterator;
    size_t queryCount;
    size_t totalMatchCount;
    size_t multipleMatchCount;
    size_t perfectMatchCount;
    size_t selectedMatchCount;

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

    vector<size_t> closestKmers;
    vector<QueryInfo> queryInfos;
    unordered_map<TaxID, unsigned int> taxCounts;


    const static uint64_t MARKER = ~0 & ~16777215;
    uint8_t hammingLookup[8][8]= {
            {0, 1, 1, 1, 2, 1, 3, 3},
            {1, 0, 1, 1, 2, 2, 3, 2},
            {1, 1, 0, 1, 2, 2, 2, 3},
            {1, 1, 1, 0, 1, 2, 3, 3},
            {2, 2, 2, 1, 0, 1, 4, 4},
            {1, 2, 2, 2, 1, 0, 4, 4},
            {3, 3, 2, 3, 4, 4, 0, 1},
            {3, 2, 3, 3, 4, 4, 1, 0}}; /// 4 means that there is no case where that value is used.


    uint8_t getHammingDistance(uint64_t kmer1, uint64_t kmer2);
    void linearSearchParallel(QueryKmer * queryKmerList, size_t & queryKmerCnt, const MmapedData<uint16_t> & targetDiffIdxList,
                              const MmapedData<TargetKmerInfo> & targetInfoList, const MmapedData<DiffIdxSplit> & diffIdxSplits,
                              Buffer<Match> & matchBuffer, const vector<int> & taxIdList, const vector<int> & speciesTaxIdList, const vector<TaxID> & genusTaxIdList,
                              FILE * matchFile);
    TaxID match2LCA(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
                    size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent);
    TaxID match2LCA2(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
                    size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent);

    static bool compareForLinearSearch(const QueryKmer & a, const QueryKmer & b);
    static bool compareConsecutiveMatches(const ConsecutiveMatches & a, const ConsecutiveMatches & b);
    static bool compareConsecutiveMatches2(const ConsecutiveMatches & a, const ConsecutiveMatches & b);
    static bool compareConsecutiveMatches3(const ConsecutiveMatches & a, const ConsecutiveMatches & b);
    void fillQueryKmerBufferParallel(QueryKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt, Query * queryList);
    TaxID chooseBestTaxon(NcbiTaxonomy & ncbiTaxonomy, const size_t & queryLength, const int & currentQuery, const size_t & offset, const size_t & end, Match * matchList, Query * queryList);

    void writeReadClassification(Query * queryList, int queryNum , ofstream & readClassificationFile);
    void writeReportFile(const char * queryFileName, NcbiTaxonomy & ncbiTaxonomy, const int numOfQuery);
    void analyseResultParallel(NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments, char * matchFileName, int seqNum, Query * queryList);
    void writeReport(FILE * fp, const NcbiTaxonomy & ncbiTaxonomy, const unordered_map<TaxID, TaxonCounts> & cladeCounts, unsigned long totalReads,TaxID taxID = 0, int depth = 0);
    unsigned int cladeCountVal(const std::unordered_map<TaxID, TaxonCounts>& map, TaxID key);
    void compareDna(uint64_t & query, vector<uint64_t> & targetKmersToCompare, const size_t & startIdx, vector<size_t> & selectedMatches, vector<uint8_t> & selectedHamming);
    void writeMatches(Buffer<Match> & matchBuffer, FILE * matchFile);
    static bool compareForWritingMatches(const Match & a, const Match & b);
    static bool sortByTaxId(const Match & a, const Match & b);
    void findConsecutiveMatches(vector<ConsecutiveMatches> & list, Match * matchList, size_t end, size_t begin);
    void getBestGenusLevelMatchCombination(vector<ConsecutiveMatches> & chosenMatchCombination, Match * matchList, size_t end, size_t offset);

    void getMatchCombinationForCurGenus(vector<ConsecutiveMatches> & coMatches, vector<vector<ConsecutiveMatches>> & genus, Match * matchList);
    void getMatchCombinationForCurGenus2(vector<ConsecutiveMatches> & coMatches, vector<vector<ConsecutiveMatches>> & genus, Match * matchList);


    void getTheBestGenus(vector<vector<ConsecutiveMatches>> & genus, vector<ConsecutiveMatches> & choosed);
    void getSubsets(vector<int> & subset, vector<vector<int>> & uniqueSubset, int k, int n);
    float scoreSubset(vector<ConsecutiveMatches> & subset);

    //void getQueryKmerSplits(int threadNum, )
public:
    void startClassify(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName, const char * diffIdxSplitFileName, vector<int> & taxIdList, const LocalParameters & par);
    static uint64_t getNextTargetKmer(uint64_t lookingTarget, const uint16_t * targetDiffIdxList, size_t & diffIdxPos);
    int getNumOfSplits() const;
    Classifier();
    ~Classifier();
    void performanceTest(NcbiTaxonomy & ncbiTaxonomy, Query * queryList, int numOfquery);
    void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy);
};

inline uint8_t Classifier::getHammingDistance(uint64_t kmer1, uint64_t kmer2) {
    uint8_t hammingDist = 0;
    for(int i = 0; i < 8 ; i++){
        hammingDist += hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)];
        kmer1 >>= 3U;
        kmer2 >>= 3U;
    }
    return hammingDist;
}

inline uint64_t Classifier::getNextTargetKmer(uint64_t lookingTarget, const uint16_t* targetDiffIdxList, size_t & diffIdxPos){
    uint16_t fragment;
    uint16_t check = (0x1u << 15u);
    uint64_t diffIn64bit = 0;

    fragment = targetDiffIdxList[diffIdxPos];
    diffIdxPos++;
    while (!(fragment & check)){
        diffIn64bit |= fragment;
        diffIn64bit <<= 15u;
        fragment = targetDiffIdxList[diffIdxPos];
        diffIdxPos++;
    }
    fragment &= ~check;
    diffIn64bit |= fragment;

    return diffIn64bit + lookingTarget;
}

#endif //ADKMER4_SEARCHER_H
