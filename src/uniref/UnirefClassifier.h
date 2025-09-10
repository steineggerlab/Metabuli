#ifndef METABULI_UNIREF_CLASSIFIER_H
#define METABULI_UNIREF_CLASSIFIER_H

#include <iostream>
#include <fstream>
#include "LocalParameters.h"
#include "UnirefTree.h"
#include "common.h"
#include "KmerExtractor.h"  
#include "KmerMatcher.h"
#include "Match.h"

class UnirefClassifier {
private:
    const LocalParameters & par;
    std::string dbDir;
    std::string inputFileName;
    std::string outDir;

    GeneticCode * geneticCode = nullptr;
    UnirefTree * unirefTree = nullptr;
    KmerExtractor * kmerExtractor = nullptr;
    KmerMatcher * kmerMatcher = nullptr;
    int matchPerKmer = 4;
    

    // Results
    std::unordered_map<uint32_t, size_t> clusterCnts;
    int writeCnt = 0;
    ofstream outFile;


    void analyzeMatches(
        const Buffer<Match_AA> & matchBuffer,
        std::vector<ProteinQuery> & queryList,
        size_t seqCnt);

    uint32_t assignUniref(
        const Match_AA * matchList,
        MatchBlock & block,
        std::vector<ProteinQuery> &queryList);

    void writeResults(
        const std::vector<ProteinQuery> & queryList,
        size_t seqCnt);

    size_t calculateBufferSize();

        
public:
    UnirefClassifier(const LocalParameters &par, UnirefTree * unirefTree);
    ~UnirefClassifier() {
        delete kmerExtractor;
        delete kmerMatcher;
        outFile.close();
    }
    int classify();
};
#endif //METABULI_UNIREF_CLASSIFIER_H