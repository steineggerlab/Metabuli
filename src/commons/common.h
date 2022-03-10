#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include <utility>


#include "NcbiTaxonomy.h"
#define kmerBufSize 10000000000 // 10000000000 | 286000000
#define ThreadNum 32
#define SplitNum 4096
#define PRINT false

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
    bool isClassified;
    bool newSpecies;
    int queryLength;
    int queryLength2;
    int kmerCnt;
    string name;
    unordered_map<TaxID,int> taxCnt;

    bool operator==(int id) const { return queryId == id;}
//    Query(int id, int classification_, float score, bool isClassified_, bool newSpecies, uint32_t len, string name_)
//    :queryId(id), classification(classification_), score(score), isClassified(isClassified_), newSpecies(newSpecies), queryLength(len),
//    name(std::move(name_)) { }
    Query():queryId(0), classification(0), score(0.0f), isClassified(false), queryLength(0) {}
};



#endif //ADCLASSIFIER2_COMMON_H
