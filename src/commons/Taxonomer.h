#ifndef METABULI_TAXONOMER_H
#define METABULI_TAXONOMER_H

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>

#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include "Match.h"
#include "common.h"
#include "BitManipulateMacros.h"
#include "MetamerPattern.h"
#include "EvalueComputation.h"

using namespace std;

struct TaxonScore {
    TaxID taxId;
    MatchScore score;
    int hammingDist;
    bool LCA;
    TaxonScore(TaxID taxId, MatchScore score, int hammingDist, bool LCA) :
            taxId(taxId), score(score), hammingDist(hammingDist), LCA(LCA) {}
    TaxonScore() : taxId(0), score(), hammingDist(0), LCA(false) {}
};

class Taxonomer {
private:
    const LocalParameters & par;
    TaxonomyWrapper * taxonomy;
    const MetamerPattern *metamerPattern = nullptr;
    SubstitutionMatrix * substitutionMatrix = nullptr;
    EvalueComputation *evaluer = nullptr;

    // spaced k-mer
    int unmaskedPos[9];
    int spaceNum;
    int kmerLen;
    int windowSize;
    uint32_t windowMask;

    // Parameters from user
    int maxGap;
    int accessionLevel;
    int minSSMatch;
    size_t minConsCnt;
    size_t minConsCntEuk;
    int eukaryotaTaxId;
    float tieRatio;

    // Internal
    int denominator;
    int maxCodonShift;
    int dnaShift;
    // int smerLength;
    int minSubSpeciesMatch;
    size_t dbSize;
    double logMaxEValue;
    bool useEvalueFilter = false;

    // vector<const Match *> speciesMatches;

    // chooseBestTaxon
    unordered_map<TaxID, unsigned int> taxCnt;

    // getBestSpeciesMatches
    vector<MatchPath> matchPaths;
    vector<MatchPath> combinedMatchPaths;
    vector<TaxID> maxSpecies;

    // getMatchPaths
    vector<bool> connectedToNext;
    vector<MatchPath> localMatchPaths;

    // lowerRankClassification
    unordered_map<TaxID, TaxonCounts> cladeCnt;

    // filterRedundantMatches
    const Match **bestMatchForQuotient;
    TaxID *bestMatchTaxIdForQuotient;
    uint8_t *minHammingForQuotient;
    size_t arraySize_filterRedundantMatches;


    // Output
    unordered_map<TaxID, unsigned int> taxCounts;

    void ensureArraySize(size_t newSize);

    void printSpeciesMatches (
       const Match *matchList,
       const std::pair<size_t, size_t> & bestSpeciesRange
    );

    TaxonScore getBestSpeciesMatches(
        std::pair<size_t, size_t> & bestSpeciesRange,
        const Match *matchList,
        size_t end,
        size_t offset,
        Query & query);

    void getMatchPaths(
        const Match * matchList,
        size_t matchNum,
        vector<MatchPath> & matchPaths,
        TaxID speciesId);

    void getMatchPaths2(
        const Match * matchList,
        size_t matchNum,
        vector<MatchPath> & matchPaths,
        TaxID speciesId); 

    void getMatchPaths3(
        const Match * matchList,
        size_t matchNum,
        vector<MatchPath> & matchPaths,
        TaxID speciesId); 

    MatchPath makeMatchPath(
        const Match * match
    );

    void makeMatchPath(
        const Match * match,
        size_t index
    );

    MatchScore combineMatchPaths(
        vector<MatchPath> & matchPaths,
        size_t matchPathStart,
        vector<MatchPath> & combMatchPaths,
        size_t combMatchPathStart,
        int queryLength);
        
    bool isMatchPathOverlapped(const MatchPath & matchPath1, const MatchPath & matchPath2);
    void trimMatchPath(MatchPath & path1, const MatchPath & path2, int overlapLength);
    void trimMatchPath2(MatchPath & path1, const MatchPath & path2, int overlapLength);
    void sortMatchPath(std::vector<MatchPath> & matchPaths, size_t i);

public:
    Taxonomer(
        const LocalParameters & par, 
        TaxonomyWrapper * taxonomy, 
        const MetamerPattern *metamerPattern);

    ~Taxonomer();

    void chooseBestTaxon(uint32_t currentQuery,
                         size_t offset,
                         size_t end,
                         const Match *matchList,
                         vector<Query> & queryList,
                         const LocalParameters &par);      

    void filterRedundantMatches(const Match *matchList,
                                const std::pair<size_t, size_t> & bestSpeciesRange,
                                unordered_map<TaxID, unsigned int> & taxCnt,
                                int queryLength);

    TaxID lowerRankClassification(const unordered_map<TaxID, unsigned int> & taxCnt, TaxID speciesID, int queryLength);

    void getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> & taxCnt,
                               unordered_map<TaxID, TaxonCounts> & cladeCnt,
                               TaxID spciesID);

    TaxID BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root, unsigned int maxCnt);

    // Getters
    unordered_map<TaxID, unsigned int> & getTaxCounts() { return taxCounts; }

};


#endif //METABULI_TAXONOMER_H