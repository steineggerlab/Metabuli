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
    int minConsCnt;
    int minConsCntEuk;
    int eukaryotaTaxId;
    float tieRatio;

    // Internal
    int denominator;
    vector<const Match *> speciesMatches;

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

    // lowerRankClassification
    unordered_map<TaxID, TaxonCounts> cladeCnt;

    // Output
    unordered_map<TaxID, unsigned int> taxCounts;

    depthScore DFS(const Match *curMatch,
                   const vector<const Match *> &linkedMatchesKeys,
                   const vector<const Match *> &linkedMatchesValues,
                   const vector<size_t> &linkedMatchesIndices,
                   size_t depth, size_t MIN_DEPTH,
                   unordered_map<const Match *, depthScore> &match2depthScore,
                   float score, int hammingDist);


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
    
    float combineMatchPaths(vector<MatchPath> & matchPaths,
                            size_t matchPathStart,
                            vector<MatchPath> & combMatchPaths,
                            size_t combMatchPathStart,
                            int readLength);

    bool isMatchPathOverlapped(const MatchPath & matchPath1, const MatchPath & matchPath2);

    bool isMatchPathLinked(const MatchPath & matchPath1, const MatchPath & matchPath2);

    void mergeMatchPaths(const MatchPath & source, MatchPath & target);

    void trimMatchPath(MatchPath & path1, const MatchPath & path2, int overlapLength);

    void filterRedundantMatches(vector<const Match *> & speciesMatches,
                                unordered_map<TaxID, unsigned int> & taxCnt);

    static bool isConsecutive(const Match * match1, const Match * match2);

    static bool isConsecutive_diffFrame(const Match * match1, const Match * match2);

    TaxonScore getBestSpeciesMatches(vector<const Match *> &speciesMatches, const Match *matchList, size_t end,
                                     size_t offset, int queryLength);
    
    // TaxonScore getBestSpeciesMatches(vector<Match> &speciesMatches, const Match *matchList, size_t end,
    //                                  size_t offset, Query & currentQuery);
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


    bool compareMatchPaths(const MatchPath& a, const MatchPath& b) const {
        if (a.score != b.score)
            return a.score < b.score;
        if (a.hammingDist != b.hammingDist)
            return a.hammingDist < b.hammingDist;
        return a.start < b.start;
    }
};


#endif //METABULI_TAXONOMER_H