//
// Created by KJB on 11/09/2020.
//

#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include "NcbiTaxonomy.h"
#define kmerBufSize 10000000000 // 10000000000 | 2860000000
#define ThreadNum 64
#define SplitNum 4096
#define PRINT false;

//#define kmerBufSize 1000000000
#define AApart(x) x & ()
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
    uint32_t queryLength;
    string name;
    unordered_map<TaxID,int> taxCnt;

    bool operator==(int id) const { return queryId == id;}
    Query(int id, bool isClassified_, const string & name_, int classification_, float score, uint32_t len)
    :queryId(id), isClassified(isClassified_), name(name_), classification(classification_), score(score), queryLength(len) { }
    Query():queryId(0), isClassified(false), classification(0), score(0.0f), queryLength(0) {}
};



#endif //ADCLASSIFIER2_COMMON_H
