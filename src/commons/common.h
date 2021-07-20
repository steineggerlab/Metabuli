//
// Created by KJB on 11/09/2020.
//

#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include "NcbiTaxonomy.h"
#define kmerBufSize 10000000000 // 10,000,000,000
#define ThreadNum 32
#define SplitNum 4096
//#define kmerBufSize 1000000000
#define AApart(x) x & ()
struct Sequence{
    Sequence(size_t start, size_t end, size_t length) : start(start), end(end), length(length) { }
    Sequence() : start(0), end(0), length(0) { }
    size_t start;
    size_t end;
    size_t length;
};

struct Query{
    int queryId;
    bool isClassified; //필요한가?
    string name;
    int classification;
    float coverage;
    uint32_t queryLength;
    unordered_map<TaxID,int> taxCnt;
    bool operator==(int id) const { return queryId == id;}
    Query(int id, bool isClassified_, const string & name_, int classification_, float coverage_, uint32_t len)
    :queryId(id), isClassified(isClassified_), name(name_), classification(classification_), coverage(coverage_), queryLength(len) { }
    Query():queryId(0), isClassified(false), classification(0), coverage(0.0f), queryLength(0) {}
};



#endif //ADCLASSIFIER2_COMMON_H
