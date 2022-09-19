#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include <utility>
#include "NcbiTaxonomy.h"
#include <iostream>
#define kmerBufSize 286'000'000 // 10'000'000'000  //1'000'000'000 //  10'000'000'000 | 286'000'000 (16 byte x 1 giga = 16 GB)
                                // 10'000'000'000 -> build_dir 397G RAM
#define SplitNum 4096
#define PRINT true

#define likely(x) __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)

struct Sequence{
    Sequence(size_t start, size_t end, size_t length) : start(start), end(end), length(length) { }
    Sequence() : start(0), end(0), length(0) { }
    size_t start;
    size_t end;
    int length;
};

struct Query{
    int queryId;
    int classification;
    float score;
    int queryLength;
    int queryLength2;
    int kmerCnt;
    bool isClassified;
    bool newSpecies;

    std::string name;
    std::unordered_map<TaxID,int> taxCnt;

    bool operator==(int id) const { return queryId == id;}
//    Query(int id, int classification_, float score, bool isClassified_, bool newSpecies, uint32_t len, string name_)
//    :queryId(id), classification(classification_), score(score), isClassified(isClassified_), newSpecies(newSpecies), queryLength(len),
//    name(std::move(name_)) { }
    Query():queryId(0), classification(0), score(0.0f), queryLength(0), queryLength2(0), isClassified(false), newSpecies(false) {}
};



#endif //ADCLASSIFIER2_COMMON_H
