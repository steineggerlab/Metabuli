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

#include <vector>
#include <algorithm>
#include "FastSort.h"
#include "kseq.h"
#include "KSeqBufferReader.h"



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

    int numOfSplit;
    SeqIterator * seqIterator;
    size_t queryCount;
    size_t totalMatchCount;
    size_t multipleMatchCount;
    size_t perfectMatchCount;
    size_t closestCount;
//    uint64_t currentTargetKmer;
//    uint64_t currentQuery;
//    uint64_t nextTargetKmer;
    size_t currentTargetPos;
    int isMatched;
    vector<MatchedKmer> matchedKmerList;
    vector<size_t> closestKmers;
    int currentHamming;
    ExtractStartPoint ESP;


    const static uint64_t MARKER = ~0 & ~16777215;
    uint8_t hammingLookup[8][8]= {
            {0, 1, 1, 1, 2, 1, 3, 3},
            {1, 0, 1, 1, 2, 2, 3, 2},
            {1, 1, 0, 1, 2, 2, 2, 3},
            {1, 1, 1, 0, 1, 2, 3, 3},
            {2, 2, 2, 1, 0, 1, 4, 4},
            {1, 2, 2, 2, 1, 0, 4, 4},
            {3, 3, 2, 3, 4, 4, 0, 1},
            {3, 2, 3, 3, 4, 4, 1, 0}};
    uint8_t getHammingDistance(uint64_t kmer1, uint64_t kmer2);
    static uint64_t getNextTargetKmer(uint64_t lookingTarget, const uint16_t * targetDiffIdxList, size_t & diffIdxPos);
    void linearSearch(Kmer * queryKmerList, size_t & numOfKmer, const MmapedData<uint16_t> & targetDiffIdxList, const MmapedData<KmerInfo> & targetInfoList, const vector<int> & taxIdList);
    void writeResultFile(vector<MatchedKmer> & matchList, const char * queryFileName);
    static TaxID selectALeaf(unordered_map<TaxID, int> & listOfLCAs, NcbiTaxonomy & ncbiTaxonomy, const size_t & length, float coverageThr = 0.8);
    TaxID selectLcaFromTaxIdList(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
                           size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent);
//    TaxID selectTaxForSet(const std::vector<int> &setTaxa, NcbiTaxonomy const *taxonomy, const float majorityCutoff,
//                          size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent);


public:
    void startClassify(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName, const vector<int> & taxIdList);
    void startClassify2(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName, const vector<int> & taxIdList);
    void analyseResult(const char * queryFileName, NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments);
    void analyseResult2(NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments);
    int getNumOfSplits() const;
    Classifier();
    ~Classifier();
};
#endif //ADKMER4_SEARCHER_H
