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
#include <cstdio>
#include <ctime>
#include <vector>
#include <algorithm>
#include <regex>
#include "FastSort.h"
#include "KSeqWrapper.h"
#include "LocalParameters.h"
#include <set>
#include <cmath>
#include "Match.h"
#include <unordered_set>
#include "LocalUtil.h"
#include "QueryIndexer.h"
#include "ReducedKmerMatcher.h"
#include "KmerExtractor.h"
#include "Taxonomer.h"
#include "Reporter.h"
#include <cassert>
#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;

// new code
// void saveQueryIdToFile(const std::vector<Query>& queryIdMap, const std::string& queryIdFileDir);
// void loadQueryIdFromFile(const std::string& queryIdFileDir, 
//                            std::unordered_map<std::string, std::string>& queryIdMap);
void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t &localBufIdx);
void getDiffIdx(const uint64_t &lastKmer, const uint64_t &entryToWrite, FILE *handleKmerTable, uint16_t *buffer, size_t &localBufIdx);
void writeDiffIdx(uint16_t *buffer, FILE *handleKmerTable, uint16_t *toWrite, size_t size, size_t &localBufIdx);
void writeQueryKmerFile(QueryKmerBuffer * queryKmerBuffer, const std::string& queryKmerFileDir);

void processKmerQuery(const std::string& queryKmerFileDir, 
                      const std::string& groupFileDir,
                      size_t processedReadCnt);

void makeGraph(const std::string& queryKmerFileDir, 
               std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> & relation);
               
void makeGroups(const std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> & relation, //std::unordered_map<std::string, std::unordered_map<std::string, int>>& relation,
                std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo,
                std::vector<uint32_t> &queryGroupInfo);

void saveGroupsToFile(const std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo, 
                      const std::vector<uint32_t>& queryGroupInfo,
                      const std::string& groupFileDir);

void loadGroupInfo(const std::string& groupFileDir, 
                    std::unordered_map<std::string, std::unordered_set<std::string>>& groupInfo);

void loadQueryGroupInfo(const std::string& groupFileDir, 
                        std::vector<uint32_t> & queryGroupInfo);

                        

void loadMetabuliResult(const std::string& filename, 
                        std::unordered_map<std::string, std::string>& queryIdMap, 
                        std::unordered_map<std::string, int>& metabuliResult);

void getRepLabel(const std::unordered_map<std::string, int>& metabuliResult,
                 const std::unordered_map<std::string, std::unordered_set<std::string>>& groupInfo,
                 std::unordered_map<std::string, int>& repLabel);

void applyRepLabel(const std::string& resultFileDir, 
                   const std::string& newResultFileDir, 
                   const std::vector<uint32_t>& queryGroupInfo,
                   const std::unordered_map<uint32_t, int>& repLabel);

extern int kmerQueryFileNumber;
const size_t bufferSize = 1024 * 1024;

// DisjointSet class
class DisjointSet {
public:
    // std::unordered_map<std::string, std::string> parent;
    // std::unordered_map<std::string, int> rank;

    std::unordered_map<uint32_t, uint32_t> parent;
    std::unordered_map<uint32_t, int> rank;

    void makeSet(uint32_t item);

    uint32_t find(uint32_t item);

    void unionSets(uint32_t set1, uint32_t set2);
};

class Classifier {
protected:
    // Parameters
    string dbDir;
    size_t matchPerKmer;

    // Agents
    QueryIndexer * queryIndexer;
    KmerExtractor * kmerExtractor;
    KmerMatcher * kmerMatcher;
    Taxonomer * taxonomer;
    Reporter * reporter;
    NcbiTaxonomy * taxonomy;

public:
    void startClassify(const LocalParameters &par);

    explicit Classifier(LocalParameters & par);

    virtual ~Classifier();

};


#endif //ADKMER4_SEARCHER_H
