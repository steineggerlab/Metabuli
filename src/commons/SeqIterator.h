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
#include "KmerBuffer.h"
#include <algorithm>
#include "omp.h"
#include "kseq.h"
#include "KSeqBufferReader.h"
#include "ProdigalWrapper.h"

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define kmerLength 8

#define nuc2int(x) (x & 6u)>>1u

using namespace std;

typedef struct PredictedGene {
    PredictedGene(int start, int end) : start(start), end(end) { }
    int start;
    int end;
}PredictedGene;

typedef struct ExtractStartPoint {
    uint32_t frame;
    uint32_t startOfFrame;
}ExtractStartPoint;

class SeqIterator
{
private:
    string iRCT;
    vector<int> aaFrames[6];
    uint64_t powers[10];
    int nuc2aa[4][4][4];
    uint64_t nuc2num[4][4][4];

    uint64_t addDNAInfo(uint64_t & kmer, const string&, int forOrRev, int startOfKmer, int frame);
    void addDNAInfo2(uint64_t & kmer, Sequence & seq, MmapedData<char> & seqFile, const int & forOrRev, const int & startOfKmer, const int & frame);
    //void numAA2alphabet();

public:
    ExtractStartPoint fillKmerBuffer(const string & seq, Kmer* kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP);
    ExtractStartPoint fillKmerBuffer2(Sequence seq, MmapedData<char> & seqFile, Kmer * kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP);
    void fillKmerBuffer3(const string & seq , KmerBuffer & kmerBuffer, size_t & posToWrite, const int & seqID);
    size_t whatNameWouldBeGood(KmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt);
    size_t whatNameWouldBeGoodWithFramePrediction(KmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt);
    string reverseCompliment(string & read) const ;
    void dna2aa(const string& seq);
    void dna2aa2(const Sequence & seq, const MmapedData<char> & seqFile);
    void getSeqSegmentsWithoutHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile);
    void getSeqSegmentsWithHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile);
    size_t getNumOfKmerForSeq(const string & seq);
    SeqIterator();
};
#endif //ADKMER4_KMEREXTRACTOR_H

