
#ifndef GROUP_GENERATOR_H
#define GROUP_GENERATOR_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "IndexCreator.h"
#include "SeqIterator.h"
#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include "QueryIndexer.h"
#include "ReducedKmerMatcher.h"
#include "KmerExtractor.h"
#include "Taxonomer.h"
#include "Reporter.h"
#include "KSeqWrapper.h"
#include "DeltaIdxReader.h"
#include <set>
#include <cassert>
#include <thread>
#include <atomic>

#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;

struct Relation {
    uint32_t id1;
    uint32_t id2;
    uint32_t weight;

    bool operator==(const Relation& other) const {
        return id1 == other.id1 && id2 == other.id2;
    }
};

// struct compareForKmerFilter {
//     bool operator()(const QueryKmer& a, const QueryKmer& b) const {
//         if (a.info.sequenceID != b.info.sequenceID)
//             return a.info.sequenceID < b.info.sequenceID;
//         return a.info.pos < b.info.pos;
//     }
// };

struct relation_hash {
    size_t operator()(const Relation& r) const {
        return hash<uint64_t>()((static_cast<uint64_t>(r.id1) << 32) | r.id2);
    }
};

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

struct MetabuliInfo {
    int label;
    float score;
    string name;
};

class KmerFileHandler {
private:
    std::vector<uint64_t> kmerBoundaries;
    bool boundariesInitialized = false;
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
    vector<pair<uint64_t, QueryKmerInfo>> getNextKmersBatch(const MmapedData<uint16_t>& diffList,
                                                            const MmapedData<QueryKmerInfo>& infoList,
                                                            size_t& idxDiff,
                                                            size_t& idxInfo,
                                                            uint64_t& currentVal,
                                                            size_t maxBatchSize);

    void writeQueryKmerFile(
        Buffer<Kmer>& queryKmerBuffer, 
        const string& outDir, 
        size_t& numOfSplits, 
        size_t numOfThreads, 
        size_t processedReadCnt);

    std::vector<std::pair<size_t, size_t>> getKmerRanges(
        const Buffer<Kmer> & kmerBuffer, 
        size_t offset);

    void writeQueryKmerFile2(
        Buffer<Kmer>& queryKmerBuffer, 
        const string& outDir,
        size_t& numOfSplits, 
        size_t numOfThreads, 
        size_t processedReadCnt);

    
    
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
    const LocalParameters & par;
    string commonKmerDB;
    string taxDbDir;
    string orgRes;
    string outDir;
    size_t matchPerKmer;
    int kmerFormat;
    
    // Agents    
    GeneticCode * geneticCode;
    QueryIndexer *queryIndexer;
    KmerExtractor *kmerExtractor;
    Reporter *reporter;    
    // KmerMatcher * kmerMatcher;
    TaxonomyWrapper *taxonomy;
    KmerFileHandler *kmerFileHandler;
    
    unordered_map<TaxID, TaxID> taxId2speciesId;
    unordered_map<TaxID, TaxID> taxId2genusId;

public:
    GroupGenerator(LocalParameters & par);

    void startGroupGeneration(const LocalParameters & par);
    
    void filterCommonKmers(Buffer<Kmer>& queryKmerBuffer,
                           Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
                           const string & db="");
    void filterCommonKmers(Buffer<Kmer>& queryKmerBuffer,
                                       const string & db);

    void makeGraph(
        size_t &numOfSplits, 
        size_t &numOfThreads, 
        size_t &numOfGraph,
        size_t processedReadCnt);

    void makeGraph2(
        size_t &numOfSplits, 
        size_t &numOfThreads, 
        size_t &numOfGraph,
        size_t processedReadCnt);

    void saveSubGraphToFile(
        const unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> &subRelation, 
        const size_t counter_now);

    double mergeRelations(size_t numOfGraph,
                          const vector<MetabuliInfo>& metabuliResult,
                          const double thresholdK);

    void mergeRelations(size_t numOfGraph);

    void makeGroups(int groupKmerThr,
                    unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                    vector<int> &queryGroupInfo);

    void saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                          const vector<int> &queryGroupInfo,
                          const vector<MetabuliInfo>& metabuliResult);

    void loadGroupsFromFile(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                       vector<int> &queryGroupInfo,
                       const string &groupFileDir,
                       const string &jobId);

    void loadMetabuliResult(vector<MetabuliInfo>& metabuliResult);

    void getRepLabel(
        vector<MetabuliInfo>& metabuliResult, 
        const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
        unordered_map<uint32_t, int> &repLabel, 
        const float groupScoreThr);
    
    void loadRepLabel(std::unordered_map<uint32_t, int> &repLabel);

    void applyRepLabel( 
        const vector<int> &queryGroupInfo, 
        const unordered_map<uint32_t, int> &repLabel, 
        const float groupScoreThr);
    
    void makeGroupsFromBinning(const string &binningFileDir, 
                               unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> &relation,
                               unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                               vector<int> &queryGroupInfo, 
                               int groupKmerThr);

    void tempFunction();

    ~GroupGenerator();
};


#endif // GROUP_GENERATOR_H