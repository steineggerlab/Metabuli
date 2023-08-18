#ifndef METABULI_FILTERER_H
#define METABULI_FILTERER_H

#include "LocalUtil.h"
#include "QueryIndexer.h"
#include "ReducedKmerMatcher.h"
class QueryFilter {
private:
    QueryIndexer * queryIndexer;
    KmerMatcher * kmerMatcher;

    std::string in1, in2, out1, out2, reportFileName; // input and output file names

    void setInputAndOutputFiles(const LocalParameters & par);

public:
    void filterReads(LocalParameters & par);
    explicit QueryFilter(LocalParameters & par);
    ~QueryFilter();
};


#endif //METABULI_FILTERER_H
