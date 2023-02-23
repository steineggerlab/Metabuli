#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include <utility>
#include "NcbiTaxonomy.h"
#include <iostream>
#define kmerBufSize 1'000'000'000  // 1'000'000'000 //286'000'000 //   10'000'000'000 | 286'000'000 (16 byte x 1 giga = 16 GB)
// 10'000'000'000 -> build_dir 397G RAM
// 1'000'000'000 -> build_dir 39.7G RAM
#define SplitNum 4096
#define PRINT false
#define likely(x) __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)

struct Sequence{
    Sequence(size_t start, size_t end, size_t length, size_t seqLength = 0)
            : start(start), end(end), length(length), seqLength(seqLength) {}
    Sequence() : start(0), end(0), length(0), seqLength(0) { }
    size_t start;
    size_t end;
    size_t length;
    size_t seqLength;
};

struct Query{
    int queryId;
    int classification;
    float score;
    float coverage;
    int hammingDist;
    int queryLength;
    int queryLength2;
    int kmerCnt;
    bool isClassified;
    bool newSpecies; // 36 byte

    std::string name;
    std::unordered_map<TaxID,int> taxCnt; // 8 byte per element

    bool operator==(int id) const { return queryId == id;}

    Query(int queryId, int classification, float score, float coverage, int hammingDist, int queryLength,
          int queryLength2, int kmerCnt, bool isClassified, bool newSpecies, std::string name)
            : queryId(queryId), classification(classification), score(score), coverage(coverage),
              hammingDist(hammingDist), queryLength(queryLength), queryLength2(queryLength2), kmerCnt(kmerCnt),
              isClassified(isClassified), newSpecies(newSpecies), name(std::move(name)) {}

    Query() : queryId(0), classification(0), score(0), coverage(0), hammingDist(0), queryLength(0),
              queryLength2(0), kmerCnt(0), isClassified(false), newSpecies(false) {}
};

inline bool fileExist(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

void process_mem_usage(double& vm_usage, double& resident_set);



#endif //ADCLASSIFIER2_COMMON_H
