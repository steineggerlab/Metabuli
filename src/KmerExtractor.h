//
// Created by KJB on 01/09/2020.
//

#ifndef ADKMER4_KMEREXTRACTOR_H
#define ADKMER4_KMEREXTRACTOR_H
#include <iostream>
#include <vector>
#include "Kmer.h"
#include "printBinary.h"
#include "common.h"
#define kmerLength 8

#define nuc2int(x) (x & 6u)>>1u

using namespace std;

typedef struct ExtractStartPoint {
    size_t frame;
    size_t startOfFrame;
}ExtractStartPoint;

class KmerExtractor
{
private:
    string iupacReverseComplementTable;
    vector<int> aaFrames[6];
    uint64_t powers[10];
    int nuc2aa[4][4][4];
    uint64_t nuc2num[4][4][4];

    uint64_t addDNAInfo(uint64_t, const string&, int i2, int i3);

public:
    ExtractStartPoint fillKmerBuffer(const string * dnaSeq, Kmer* kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP);
    string reverseCompliment(string & read) const ;
    void dna2aa(const string& forward, const string & reverse);
    KmerExtractor();
};
#endif //ADKMER4_KMEREXTRACTOR_H

