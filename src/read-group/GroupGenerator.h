
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
    uint16_t weight;

    Relation() : id1(0), id2(0), weight(0) {}
    Relation(uint32_t id1, uint32_t id2, uint16_t weight) : id1(id1), id2(id2), weight(weight) {}

    bool operator==(const Relation& other) const {
        return id1 == other.id1 && id2 == other.id2;
    }

    bool operator<(const Relation& other) const {
        if (id1 != other.id1)
            return id1 < other.id1;
        if (id2 != other.id2)
            return id2 < other.id2;
        return weight < other.weight;
    }

    static bool compare(const Relation& a, const Relation& b) {
        if (a.id1 != b.id1)
            return a.id1 < b.id1;
        return a.id2 < b.id2;
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

struct OrgResult {
    int label;
    float score;
    string name;
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
    TaxonomyWrapper *taxonomy;

    // Output
    string updatedResultFileName;
    string updatedReportFileName;
    unordered_map<TaxID, unsigned int> taxCounts;
    
    unordered_map<TaxID, TaxID> taxId2speciesId;
    unordered_map<TaxID, TaxID> taxId2genusId;
    size_t numOfSplits = 0;
    size_t numOfGraph = 0;
    std::vector<uint64_t> kmerBoundaries;
    bool boundariesInitialized = false;
    bool useOnlyTrueRelations = false; // for debug

public:
    GroupGenerator(LocalParameters & par);

    void startGroupGeneration(const LocalParameters & par);
    
    void filterCommonKmers(Buffer<Kmer>& queryKmerBuffer,
                           Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
                           const string & db="");

    void filterCommonKmers2(Buffer<Kmer>& queryKmerBuffer,
                           Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
                           const string & db="");

    void makeGraph(size_t processedReadCnt);
    
    void saveSubGraphToFile(
        const unordered_map<uint64_t, uint16_t> & pair2weight,
        const size_t counter_now);
    
    void mergeTrueRelations(
        const vector<OrgResult>& metabuliResult);

    void mergeRelations();

    void makeGroups(int groupKmerThr,
                    unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                    vector<int> &queryGroupInfo);

    void makeGroupsFromSubGraphs(
        uint32_t groupKmerThr,
        unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
        vector<int> &queryGroupInfo,
        const vector<OrgResult>& metabuliResult);

    void saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                          const vector<int> &queryGroupInfo,
                          const vector<OrgResult>& metabuliResult);

    void loadGroupsFromFile(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                       vector<int> &queryGroupInfo,
                       const string &groupFileDir,
                       const string &jobId);

    void loadOrgResult(vector<OrgResult>& metabuliResult);

    void getRepLabel(
        vector<OrgResult>& metabuliResult, 
        const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
        unordered_map<uint32_t, int> &repLabel);
    
    void loadRepLabel(std::unordered_map<uint32_t, int> &repLabel);

    void applyRepLabel( 
        const vector<int> &queryGroupInfo, 
        const unordered_map<uint32_t, int> &repLabel);

    void writeKmers(
        Buffer<Kmer>& queryKmerBuffer, 
        size_t processedReadCnt);

    std::vector<std::pair<size_t, size_t>> getKmerRanges(
        const Buffer<Kmer> & kmerBuffer, 
        size_t offset);

    ~GroupGenerator();
};


#endif // GROUP_GENERATOR_H