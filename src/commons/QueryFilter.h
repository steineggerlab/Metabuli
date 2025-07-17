#ifndef METABULI_FILTERER_H
#define METABULI_FILTERER_H

#include "LocalUtil.h"
#include "QueryIndexer.h"
#include "ReducedKmerMatcher.h"
#include "KmerExtractor.h"
#include "Taxonomer.h"
#include "Reporter.h"
#include <cstdint>

class QueryFilter {
private:
    // Parameters
    std::string dbDir;
    size_t matchPerKmer;
    int printMode;
    int seqMode;
    std::vector<std::string> contams;

    // Agents
    GeneticCode * geneticCode;
    QueryIndexer * queryIndexer;
    KmerExtractor * kmerExtractor;
    KmerMatcher * kmerMatcher;
    Taxonomer * taxonomer;
    Reporter * reporter;
    TaxonomyWrapper * taxonomy;

    // Kseq
    KSeqWrapper* filter_kseq1;
    KSeqWrapper* filter_kseq2;

    std::string in1, in2;
    std::string f1, f2, rm1, rm2; // input and output file names
    std::string readClassificationFileName;
    std::string reportFileName;
    
    bool * isFiltered;
    size_t readCounter;
    FILE * f1_fp, * f2_fp, * rm1_fp, * rm2_fp;

    void setInputAndOutputFiles(const LocalParameters & par);

    void recordFilteredReads(const vector<Query> & queryList);
    
    void printFilteredReads();

public:
    void filterReads(LocalParameters & par);
    explicit QueryFilter(LocalParameters & par);
    ~QueryFilter();
};


#endif //METABULI_FILTERER_H
