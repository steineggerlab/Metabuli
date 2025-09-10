#ifndef METABULI_QUERYINDEXOR_H
#define METABULI_QUERYINDEXOR_H

#include "LocalParameters.h"
#include "Kmer.h"
#include "Match.h"
#include "KSeqWrapper.h"
#include "LocalUtil.h"
#include "Debug.h"
#include "common.h"
#include <unordered_map>
#include <cstdint>

struct QuerySplit {
    size_t start;
    size_t end;
    size_t kmerCnt;
    size_t readCnt;

    QuerySplit(size_t start, size_t end, size_t kmerCnt, size_t readCnt) 
        : start(start), end(end), kmerCnt(kmerCnt), readCnt(readCnt) {}
};

// Input
// 1. A set of reads

// Output
// 1. size_t numOfSeq;
// 2. vector<QuerySplit> querySplits;

class QueryIndexer {
private:
    // Input
    std::string queryPath_1;
    std::string queryPath_2;
    size_t seqMode;
    // size_t matchPerKmer;
    size_t maxRam;
    size_t threads;
    int spaceNum;

    // Internal
    size_t availableRam;
    size_t bytesPerKmer;
    int kmerLen;

    // Output
    std::size_t readNum_1;
    std::size_t readNum_2;
    std::vector<QuerySplit> querySplits;
    std::size_t totalReadLength;



public:
    explicit QueryIndexer(const LocalParameters & par);
    ~QueryIndexer() = default;

    void indexQueryFile(size_t processedNum);

    // Getters
    size_t getReadNum_1() const;
    size_t getReadNum_2() const;
    const std::vector<QuerySplit> & getQuerySplits() const;
    std::size_t getTotalReadLength() const;
    size_t getAvailableRam() const;


    // Setters
    void setKmerLen(int kmerLen2) { this->kmerLen = kmerLen2; }
    void setAvailableRam();
    void setBytesPerKmer(size_t matchPerKmer) {
        bytesPerKmer = sizeof(Kmer) + matchPerKmer * sizeof(Match);
    }

};


#endif //METABULI_QUERYINDEXOR_H
