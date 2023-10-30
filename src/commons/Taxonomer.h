#ifndef METABULI_TAXONOMER_H
#define METABULI_TAXONOMER_H
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "Match.h"
#include "common.h"
#include "BitManipulateMacros.h"
#include <unordered_set>

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

class Taxonomer {
private:
    NcbiTaxonomy * taxonomy;

    // spaced k-mer
    int unmaskedPos[9];
    int spaceNum;

    // Parameters
    int maxGap;
    int minCoveredPos;
    int accessionLevel;
    int minSSMatch;

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
    
    void chooseBestTaxon2(uint32_t currentQuery,
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

    TaxonScore getBestGenusMatches(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end,
                                   size_t offset, int queryLength, const LocalParameters &par);

    TaxonScore getBestGenusMatches(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end, size_t offset,
                                   int readLength1, int readLength2, const LocalParameters &par);

    TaxonScore getBestSpeciesMatches(vector<Match> &matchesForMajorityLCA, const Match *matchList, size_t end,
                                     size_t offset, int queryLength, const LocalParameters &par);

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

    // Getters
    unordered_map<TaxID, unsigned int> & getTaxCounts() { return taxCounts; }
};


#endif //METABULI_TAXONOMER_H
