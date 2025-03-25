
#ifndef GROUP_GENERATOR_H
#define GROUP_GENERATOR_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "SeqIterator.h"
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "QueryIndexer.h"
#include "ReducedKmerMatcher.h"
#include "KmerExtractor.h"
#include "Taxonomer.h"
#include "Reporter.h"
#include "KSeqWrapper.h"
#include <set>
#include <cassert>
#include <thread>
#include <atomic>

#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;

// DisjointSet class for handling union-find operations
class DisjointSet {
public:
    DisjointSet() {}

    void makeSet(uint32_t element) {
        if (parent.find(element) == parent.end()) {
            parent[element] = element;
            rank[element] = 0;
        }
    }

    uint32_t find(uint32_t element) {
        if (parent[element] != element) {
            parent[element] = find(parent[element]);
        }
        return parent[element];
    }

    void unionSets(uint32_t elem1, uint32_t elem2) {
        uint32_t root1 = find(elem1);
        uint32_t root2 = find(elem2);

        if (root1 != root2) {
            if (rank[root1] < rank[root2]) {
                parent[root1] = root2;
            } else if (rank[root1] > rank[root2]) {
                parent[root2] = root1;
            } else {
                parent[root2] = root1;
                rank[root1]++;
            }
        }
    }

    unordered_map<uint32_t, uint32_t> parent;
    unordered_map<uint32_t, uint32_t> rank;
};

class KmerFileHandler {
protected:

public:
    KmerFileHandler() {}
    static void flushKmerBuf(uint16_t *buffer, 
                             FILE *handleKmerTable, 
                             size_t &localBufIdx);
    static void getDiffIdx(const uint64_t &lastKmer, 
                           const uint64_t &entryToWrite, 
                           FILE *handleKmerTable, 
                           uint16_t *buffer, 
                           size_t &localBufIdx,
                           size_t bufferSize);
    static void writeDiffIdx(uint16_t *buffer, 
                             FILE *handleKmerTable, 
                             uint16_t *toWrite, 
                             size_t size, 
                             size_t &localBufIdx,
                             size_t bufferSize);
    static uint64_t getNextKmer(uint64_t lookingTarget, 
                                const struct MmapedData<uint16_t> & diffList, 
                                size_t & idx);
    static void writeQueryKmerFile(QueryKmerBuffer *queryKmerBuffer, 
                                   const string& queryKmerFileDir, 
                                   size_t& numOfSplits, 
                                   size_t numOfThreads, 
                                   size_t processedReadCnt, 
                                   const string & jobId);
    static size_t bufferSize;

    template <typename T>
    static size_t  loadBuffer(FILE *fp, T *buffer, size_t size) {
        return fread(buffer, sizeof(T), size, fp);
    }

    size_t  fillKmerInfoBuffer(size_t bufferSize,
                            FILE *kmerInfoFp,
                            TaxID *infoBuffer) {
        return loadBuffer(kmerInfoFp, infoBuffer, bufferSize);
    }
};



class GroupGenerator {
protected:
    string dbDir;
    size_t matchPerKmer;
    
    // Agents
    QueryIndexer *queryIndexer;
    KmerExtractor *kmerExtractor;
    KmerMatcher *kmerMatcher;
    Taxonomer *taxonomer;
    Reporter *reporter;
    TaxonomyWrapper *taxonomy;
    KmerFileHandler *kmerFileHandler;
    SeqIterator *seqIterator;
    
    unordered_map<TaxID, TaxID> taxId2speciesId;
    unordered_map<TaxID, TaxID> taxId2genusId;

public:
    GroupGenerator(LocalParameters & par);

    void startGroupGeneration(const LocalParameters & par);

    void makeGraph(const string &queryKmerFileDir, 
                   size_t &numOfSplits, 
                   size_t &numOfThreads, 
                   size_t &numOfGraph,
                   const string &jobId);

    void saveSubGraphToFile(const map<uint32_t, map<uint32_t, uint32_t>> &subRelation, 
                            const string &subGraphFileDir, 
                            const size_t counter_now,
                            const string &jobId);

    void makeGroups(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                    const string &subGraphFileDir, 
                    vector<int> &queryGroupInfo, 
                    int groupKmerThr, 
                    size_t &numOfGraph,
                    const string &jobId);

    void saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                          const vector<int> &queryGroupInfo, 
                          const string &groupFileDir, 
                          const string &jobId);

    void loadMetabuliResult(const string &resultFileDir, 
                            vector<pair<int, float>> &metabuliResult);

    void getRepLabel(const string &groupRepFileDir, 
                     vector<pair<int, float>> &metabuliResult, 
                     const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                     unordered_map<uint32_t, int> &repLabel, 
                     const string &jobId, 
                     int voteMode, 
                     float majorityThr,
                     const float groupScoreThr);

    void applyRepLabel(const string &resultFileDir, 
                       const string &newResultFileDir, 
                       const vector<int> &queryGroupInfo, 
                       const unordered_map<uint32_t, int> &repLabel, 
                       const float groupScoreThr, 
                       const string &jobId);
    
    void makeGroupsFromBinning(const string &binningFileDir, 
                               unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> &relation,
                               unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                               vector<int> &queryGroupInfo, 
                               int groupKmerThr);

    void tempFunction();

    ~GroupGenerator();
};

#endif // GROUP_GENERATOR_H