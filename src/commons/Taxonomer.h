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
    MatchPath(int start, int end, float score, int hammingDist, const Match * startMatch, const Match * endMatch) :
         start(start), end(end), score(score), hammingDist(hammingDist), startMatch(startMatch), endMatch(endMatch) {}
    MatchPath() : start(0), end(0), score(0.f), hammingDist(0), startMatch(nullptr), endMatch(nullptr) {}
    int start;
    int end;
    float score;
    int hammingDist;
    const Match * startMatch;
    const Match * endMatch;
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
    int minConsCnt;
    int minConsCntEuk;
    int eukaryotaTaxId;
    float tieRatio;

    // Internal
    int denominator;

    struct MatchBlock {
        MatchBlock(size_t start, size_t end, int id) : start(start), end(end), id(id) {}
        MatchBlock() : start(0), end(0), id(0) {}
        size_t start;
        size_t end;
        uint32_t id;
    };



    // Output
    unordered_map<TaxID, unsigned int> taxCounts;


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

    void remainConsecutiveMatches(const vector<const Match *> & curFrameMatches,
                                  vector<MatchPath> & matchPaths,
                                  TaxID genusId);
    
    float combineMatchPaths(vector<MatchPath> & matchPaths,
                           vector<MatchPath> & combinedMatchPaths,
                           int readLength);

    bool isMatchPathOverlapped(const MatchPath & matchPath1, const MatchPath & matchPath2);

    bool isMatchPathLinked(const MatchPath & matchPath1, const MatchPath & matchPath2);

    void mergeMatchPaths(const MatchPath & source, MatchPath & target);

    void trimMatchPath(MatchPath & path1, const MatchPath & path2, int overlapLength);

    void filterRedundantMatches(vector<Match> & speciesMatches,
                                unordered_map<TaxID, unsigned int> & taxCnt);

    depthScore DFS(const vector<const Match *> &matches, const Match * curMatchIdx,
                   const map<const Match *, vector<const Match *>> &linkedMatches,
                   size_t depth, size_t MIN_DEPTH, unordered_set<const Match *> &used,
                   unordered_map<const Match *, depthScore> &match2depthScore,
                   float score, int hammingDist);

    static bool isConsecutive(const Match * match1, const Match * match2);

    static bool isConsecutive_diffFrame(const Match * match1, const Match * match2);

    TaxonScore getBestGenusMatches(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end,
                                   size_t offset, int queryLength);

    TaxonScore getBestGenusMatches(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end, size_t offset,
                                   int readLength1, int readLength2);

    TaxonScore getBestSpeciesMatches(vector<Match> &speciesMatches, const Match *matchList, size_t end,
                                     size_t offset, int queryLength);
    
    TaxonScore getBestSpeciesMatches(vector<Match> &speciesMatches, const Match *matchList, size_t end,
                                     size_t offset, int readLength1, int readLength2);

    // TaxonScore getBestGenusMatches_spaced(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end, size_t offset,
    //                                       int readLength1, int readLength2);
    // TaxonScore getBestGenusMatches_spaced(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end, size_t offset,
    //                                       int readLength1);

    TaxonScore scoreTaxon(vector<const Match *> &filteredMatches,
                          TaxID taxId,
                          int queryLength);

    TaxonScore scoreTaxon(vector<const Match *> &filteredMatches,
                          TaxID taxId,
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

    TaxID lowerRankClassification(const unordered_map<TaxID, unsigned int> & taxCnt, TaxID speciesID, int queryLength);

    void getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> & taxCnt,
                               unordered_map<TaxID, TaxonCounts> & cladeCnt,
                               TaxID spciesID);

    TaxID BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root, unsigned int maxCnt);

    // Getters
    unordered_map<TaxID, unsigned int> & getTaxCounts() { return taxCounts; }
};


#endif //METABULI_TAXONOMER_H
