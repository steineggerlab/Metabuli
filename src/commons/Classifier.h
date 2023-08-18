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
#include "KmerBuffer.h"
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
#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;



class Classifier {
protected:
    // Parameters
    string dbDir;
    size_t matchPerKmer;

    // Agents
    QueryIndexer * queryIndexer;
    KmerExtractor * kmerExtractor;
    KmerMatcher * kmerMatcher;
    Taxonomer * taxonomer;
    Reporter * reporter;
    NcbiTaxonomy * taxonomy;




public:
    void startClassify(const LocalParameters &par);

    explicit Classifier(LocalParameters & par);

    virtual ~Classifier();


};




//inline uint64_t
//Classifier::getNextTargetKmer(uint64_t lookingTarget, const uint16_t *targetDiffIdxList, size_t &diffIdxPos) {
//    uint16_t fragment;
//    uint16_t check = (0x1u << 15u);
//    uint64_t diffIn64bit = 0;
//    fragment = targetDiffIdxList[diffIdxPos];
//    diffIdxPos++;
//    while (!(fragment & check)) { // 27 %
//        diffIn64bit |= fragment;
//        diffIn64bit <<= 15u;
//        fragment = targetDiffIdxList[diffIdxPos];
//        diffIdxPos++;
//    }
//    fragment &= ~check; // not; 8.47 %
//    diffIn64bit |= fragment; // or : 23.6%
//
//    return diffIn64bit + lookingTarget;
//}




#endif //ADKMER4_SEARCHER_H
