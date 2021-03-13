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


#define AminoAcid(x) (size_t)(x & (~0 & ~16777215))
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
        ConsecutiveMatches(uint32_t begin, uint32_t end, uint32_t hamming, uint32_t gapCnt, size_t bi, size_t ei) : begin(begin), end(end), hamming(hamming), gapCnt(gapCnt), beginIdx(bi), endIdx(ei) {}
        uint32_t begin;
        uint32_t end;
        uint32_t hamming;
        uint32_t gapCnt;
        size_t beginIdx;
        size_t endIdx;
    }ConsecutiveMatches;

    typedef struct QueryInfo{
        int queryId;
        bool isClassified;
        string name;
        int taxId;
        float coverage;
        map<TaxID,int> taxCnt; ///how about using it for REPORTFILE? --> k-mer count
        size_t queryLength;
        QueryInfo(int queryId, bool isClassified, string name, int taxId, float coverage, size_t queryLength): queryId(queryId), isClassified(isClassified), name(name), taxId(taxId), coverage(coverage), queryLength(queryLength) {}
        bool operator == (const int Id) const{
            if(Id == queryId)
                return true;
            return false;
        }
    }QueryInfo;


    int numOfSplit;
    SeqIterator * seqIterator;
    size_t queryCount;
    size_t totalMatchCount;
    size_t multipleMatchCount;
    size_t perfectMatchCount;
    size_t closestCount;

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

    vector<MatchedKmer> matchedKmerList;
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

    void linearSearch(QueryKmer * queryKmerList, size_t & numOfQuery, const MmapedData<uint16_t> & targetDiffIdxList, const MmapedData<TargetKmerInfo> & targetInfoList, const vector<int> & taxIdList, const vector<int> & taxIdListAtRank);
    void writeResultFile(vector<MatchedKmer> & matchList, const char * queryFileName);
    TaxID match2LCA(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
                    size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent);
    static bool compareForAnalyzing( const MatchedKmer & a, const MatchedKmer & b);
    static bool compareForLinearSearch(const QueryKmer & a, const QueryKmer & b);
    static bool compareConsecutiveMatches(const ConsecutiveMatches & a, const ConsecutiveMatches & b);
    void fillQueryKmerBufferParallel(QueryKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt);
    TaxID chooseBestTaxon(NcbiTaxonomy & ncbiTaxonomy, const size_t & queryLength, const int & currentQuery, const size_t & offset, const size_t & end);
    static void checkAndGive(vector<uint32_t> & posList,  vector<uint8_t> & hammingList, const uint32_t & pos, const uint8_t & hammingDist);
    void writeReadClassification(vector<QueryInfo> & queryInfos, ofstream & readClassificationFile);
    void writeReportFile(const char * queryFileName, NcbiTaxonomy & ncbiTaxonomy, const int numOfQuery);
    static uint64_t getNextTargetKmer(uint64_t lookingTarget, const uint16_t * targetDiffIdxList, size_t & diffIdxPos);
    void analyseResult(NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments);
    void writeReport(FILE * fp, const NcbiTaxonomy & ncbiTaxonomy, const unordered_map<TaxID, TaxonCounts> & cladeCounts, unsigned long totalReads,TaxID taxID = 0, int depth = 0);
    unsigned int cladeCountVal(const std::unordered_map<TaxID, TaxonCounts>& map, TaxID key);
    void storeKmerMatches();
   // void getTaxCounts(const vector<QueryInfo> & queryInfos, )
public:
    void startClassify(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName, vector<int> & taxIdList, const LocalParameters & par);
    int getNumOfSplits() const;
    Classifier();
    ~Classifier();


    void performanceTest(NcbiTaxonomy & ncbiTaxonomy);
    void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy);
};
#endif //ADKMER4_SEARCHER_H
