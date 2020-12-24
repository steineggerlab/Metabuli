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

typedef struct PredictedBlock {
    PredictedBlock(int start, int end) : start(start), end(end) { }
    int start;
    int end;
    int strand; //true for forward
}PredictedBlock;

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

    uint64_t addDNAInfo_QueryKmer(uint64_t & kmer, const char * seq, int forOrRev, const int kmerCnt, const int frame);
    void getTranslationBlocks(struct _gene * genes, struct _node * nodes, PredictedBlock * blocks, size_t numOfGene, size_t length, size_t & numOfBlocks);
    void addDNAInfo_TargetKmer(uint64_t & kmer, const char * seq, const PredictedBlock& block, const int startOfKmer);

public:
    void fillQueryKmerBuffer(const char * seq , KmerBuffer & kmerBuffer, size_t & posToWrite, const int & seqID);
    size_t fillQueryKmerBufferParallel(KmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt);
    size_t fillTargetKmerBuffer(TargetKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedTaxIdCnt, const vector<int> & startsOfTaxIDs, const vector<int> & seqCntOfTaxIDs);
    size_t fillTargetKmerBuffer2(TargetKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedTaxIdCnt, const vector<int> & startsOfTaxIDs, const vector<int> & seqCntOfTaxIDs);

    string reverseCompliment(string & read) const ;
    void sixFrameTranslateASeq(const char * seq);
    void translateBlock(const char* seq, PredictedBlock & block);
    void getSeqSegmentsWithoutHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile);
    void getSeqSegmentsWithHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile);
    size_t KmerNumOfSixFrameTranslation(const string & seq);
    size_t getNumOfKmerForBlock(const PredictedBlock & block);

    void fillBufferWithKmerFromBlock(const PredictedBlock & block, const char * seq, TargetKmerBuffer & kmerBuffer, size_t & posToWrite, const int & seqID);
    void printKmerInDNAsequence(uint64_t kmer);
    SeqIterator();
};
#endif //ADKMER4_KMEREXTRACTOR_H

