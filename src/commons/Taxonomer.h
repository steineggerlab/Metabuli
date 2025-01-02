#ifndef METABULI_TAXONOMER_H
#define METABULI_TAXONOMER_H
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "Match.h"
#include "common.h"
#include "BitManipulateMacros.h"
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>

using namespace std;

struct TaxonScore {
    TaxID taxId;
    float score;
    float coverage;
    int hammingDist;
    bool LCA;
    TaxonScore(TaxID taxId, float score, float coverage, int hammingDist, bool LCA) :
            taxId(taxId), score(score), coverage(coverage), hammingDist(hammingDist), LCA(LCA) {}
    TaxonScore() : taxId(0), score(0.0f), coverage(0.0f), hammingDist(0), LCA(false) {}
};

struct depthScore {
    depthScore(size_t depth, float score, int hammingDist, const Match * endMatch) : 
        depth(depth), score(score), hammingDist(hammingDist), endMatch(endMatch) {}
    depthScore() : depth(0), score(0.f), hammingDist(0), endMatch(nullptr) {}
    size_t depth;
    float score;
    int hammingDist;
    const Match * endMatch;
};

struct MatchPath {
    MatchPath(int start, int end, float score, int hammingDist, int depth, const Match * startMatch, const Match * endMatch) :
         start(start), end(end), score(score), hammingDist(hammingDist), depth(depth), startMatch(startMatch), endMatch(endMatch) {}
    MatchPath() : start(0), end(0), score(0.f), hammingDist(0), depth(0), startMatch(nullptr), endMatch(nullptr) {}
    MatchPath(const Match * startMatch) 
        : start(startMatch->qInfo.pos),
          end(startMatch->qInfo.pos + 23),
          score(startMatch->getScore()),
          hammingDist(startMatch->hamming),
          depth(1),
          startMatch(startMatch),
          endMatch(startMatch) {}
    
    int start;                // query coordinate
    int end;                  // query coordinate
    float score;
    int hammingDist;
    int depth;
    const Match * startMatch;
    const Match * endMatch;

    void printMatchPath() {
        std::cout << start << " " << end << " " << score << " " << hammingDist << " " << depth << std::endl;
    }
};

struct MatchBlock {
    MatchBlock(size_t start, size_t end, int id) : start(start), end(end), id(id) {}
    MatchBlock() : start(0), end(0), id(0) {}
    size_t start;
    size_t end;
    uint32_t id;
};

class Taxonomer {
private:
    NcbiTaxonomy * taxonomy;

    // spaced k-mer
    int unmaskedPos[9];
    int spaceNum;

    // Parameters from user
    int maxGap;
    int minCoveredPos;
    int accessionLevel;
    int minSSMatch;
    size_t minConsCnt;
    size_t minConsCntEuk;
    int eukaryotaTaxId;
    float tieRatio;

    // Internal
    int denominator;
    // vector<const Match *> speciesMatches;

    // chooseBestTaxon
    unordered_map<TaxID, unsigned int> taxCnt;

    // getBestSpeciesMatches
    vector<MatchPath> matchPaths;
    vector<MatchPath> combinedMatchPaths;
    vector<TaxID> maxSpecies;
    vector<TaxID> speciesList;
    vector<size_t> speciesPathIdx;
    vector<size_t> speciesCombPathIdx;
    vector<float> speciesScores;

    // remainConsecutiveMatches
    vector<const Match *> linkedMatchKeys;
    vector<const Match *> linkedMatchValues;
    vector<size_t> linkedMatchValuesIdx;
    unordered_map<const Match *, depthScore> match2depthScore;

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

    depthScore DFS(const Match *curMatch,
                   const vector<const Match *> &linkedMatchesKeys,
                   size_t linkedMatchKeysIdx,
                   const vector<const Match *> &linkedMatchesValues,
                   const vector<size_t> &linkedMatchesIndices,
                   size_t depth, size_t MIN_DEPTH,
                   unordered_map<const Match *, depthScore> &match2depthScore,
                   float score, int hammingDist);

    void ensureArraySize(size_t newSize);


public:
    Taxonomer(const LocalParameters & par, NcbiTaxonomy * taxonomy);
    ~Taxonomer();

    void assignTaxonomy(const Match *matchList,
                        size_t numOfMatches,
                        std::vector<Query> & queryList,
                        const LocalParameters &par);

    void chooseBestTaxon(uint32_t currentQuery,
                         size_t offset,
                         size_t end,
                         const Match *matchList,
                         vector<Query> & queryList,
                         const LocalParameters &par);

    void remainConsecutiveMatches(const Match * matchList,
                                  size_t start,
                                  size_t end,
                                  vector<MatchPath> & matchPaths,
                                  TaxID speciesId);

    void getMatchPaths(const Match * matchList,
                       size_t start,
                       size_t end,
                       vector<MatchPath> & matchPaths,
                       TaxID speciesId);                                    

    float combineMatchPaths(vector<MatchPath> & matchPaths,
                            size_t matchPathStart,
                            vector<MatchPath> & combMatchPaths,
                            size_t combMatchPathStart,
                            int readLength);

    bool isMatchPathOverlapped(const MatchPath & matchPath1, const MatchPath & matchPath2);

    void trimMatchPath(MatchPath & path1, const MatchPath & path2, int overlapLength);

    void filterRedundantMatches(const Match *matchList,
                                const std::pair<size_t, size_t> & bestSpeciesRange,
                                unordered_map<TaxID, unsigned int> & taxCnt,
                                int queryLength);

    static bool isConsecutive(const Match * match1, const Match * match2);

    // static bool isConsecutive_diffFrame(const Match * match1, const Match * match2);

    TaxonScore getBestSpeciesMatches(std::pair<size_t, size_t> & bestSpeciesRange,
                                     const Match *matchList, size_t end,
                                     size_t offset, int queryLength);

    TaxID lowerRankClassification(const unordered_map<TaxID, unsigned int> & taxCnt, TaxID speciesID, int queryLength);

    void getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> & taxCnt,
                               unordered_map<TaxID, TaxonCounts> & cladeCnt,
                               TaxID spciesID);

    TaxID BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root, unsigned int maxCnt);

    // Getters
    unordered_map<TaxID, unsigned int> & getTaxCounts() { return taxCounts; }


    bool compareMatchPaths(const MatchPath& a, const MatchPath& b) const {
        if (a.score != b.score)
            return a.score < b.score;
        if (a.hammingDist != b.hammingDist)
            return a.hammingDist < b.hammingDist;
        return a.start < b.start;
    }
};


#endif //METABULI_TAXONOMER_H