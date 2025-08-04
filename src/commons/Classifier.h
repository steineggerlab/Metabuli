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
#include <cstdint>
#include "TaxonomyWrapper.h"

#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;




class Classifier {
protected:
    // Parameters
    string dbDir;
    size_t matchPerKmer;
    int kmerFormat;

    // Agents
    GeneticCode * geneticCode;
    QueryIndexer * queryIndexer;
    KmerExtractor * kmerExtractor;
    KmerMatcher * kmerMatcher;
    // Taxonomer * taxonomer;
    Reporter * reporter;
    TaxonomyWrapper * taxonomy;

    unordered_map<TaxID, unsigned int> taxCounts;
    unordered_map<TaxID, unsigned int> taxCounts_EM;
    unordered_map<string, vector<pair<TaxID, float>>> mappingRes;
    vector<MappingRes> mappingResList;

    void countUniqueKmerPerSpecies(vector<uint32_t> & sp2uniqKmerCnt);

public:
    void startClassify(const LocalParameters &par);

    void assignTaxonomy(const Match *matchList,
                        size_t numOfMatches,
                        std::vector<Query> & queryList,
                        const LocalParameters &par);

    explicit Classifier(LocalParameters & par);

    virtual ~Classifier();

    unordered_map<TaxID, unsigned int> & getTaxCounts() { return taxCounts; }

    void em();

};


#endif //ADKMER4_SEARCHER_H
