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
#include "Mmap.h"
#define kmerLength 8

#define nuc2int(x) (x & 6u)>>1u

using namespace std;

typedef struct ExtractStartPoint {
    uint32_t frame;
    uint32_t startOfFrame;
}ExtractStartPoint;

class SeqAlterator
{
private:
    string iRCT;
    vector<int> aaFrames[6];
    uint64_t powers[10];
    int nuc2aa[4][4][4];
    uint64_t nuc2num[4][4][4];

    uint64_t addDNAInfo(uint64_t, const string&, int i2, int i3);
    void addDNAInfo2(uint64_t & kmer, SeqSegment & seq, MmapedData<char> & seqFile, const int & forOrRev, const int & startOfKmer, const int & frame);

public:
    ExtractStartPoint fillKmerBuffer(const string * dnaSeq, Kmer* kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP);
    ExtractStartPoint fillKmerBuffer2(SeqSegment seq, MmapedData<char> & seqFile, Kmer * kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP);
    string reverseCompliment(string & read) const ;
    void dna2aa(const string& forward, const string & reverse);
    void dna2aa2(const SeqSegment & seq, const MmapedData<char> & seqFile);
    void getSeqSegments(vector<SeqSegment> & seqSegs, MmapedData<char> seqFile);

    SeqAlterator();
};
#endif //ADKMER4_KMEREXTRACTOR_H

