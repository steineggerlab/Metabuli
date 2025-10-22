#ifndef METABULI_CLASSIFIER_H
#define METABULI_CLASSIFIER_H

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
    GeneticCode * geneticCode = nullptr;
    QueryIndexer * queryIndexer = nullptr;
    KmerExtractor * kmerExtractor = nullptr;
    KmerMatcher * kmerMatcher = nullptr;
    Reporter * reporter = nullptr;
    TaxonomyWrapper * taxonomy = nullptr;

    unordered_map<TaxID, unsigned int> taxCounts;

    size_t mappingResListSize = 0;

    // EM algorithm
    unordered_map<TaxID, unsigned int> emTaxCounts;
    unordered_map<TaxID, unsigned int> reclassifyTaxCounts;
    unordered_map<TaxID, double> taxProbs;
    MappingRes * mappingResList = nullptr;
    std::vector<Classification> emResults;
    std::unordered_set<TaxID> topSpeciesSet;

    void countUniqueKmerPerSpecies(vector<uint32_t> & sp2uniqKmerCnt);

    void loadMappings(const string & mappingResFileName);

    void loadOriginalResults(const string & classificationFileName, size_t seqNum);

public:
    void startClassify(const LocalParameters &par);

    void assignTaxonomy(const Match *matchList,
                        size_t numOfMatches,
                        std::vector<Query> & queryList,
                        const LocalParameters &par);

    explicit Classifier(LocalParameters & par);

    virtual ~Classifier();

    unordered_map<TaxID, unsigned int> & getTaxCounts() { return taxCounts; }

    void getTopSpecies(const std::vector<Query> & queries) {
        for (const auto & query : queries) {
            if (query.isClassified) {
                topSpeciesSet.insert(query.topSpeciesId);
            }
        }
    }

    void em(size_t totalQueryCnt);

    void reclassify(
        const std::vector<std::pair<size_t, size_t>> & queryRanges,
        const MappingRes * mappingResList,
        const vector<double> & sp2lengthFactor,
        size_t totalQueryCnt);

};


#endif //METABULI_CLASSIFIER_H
