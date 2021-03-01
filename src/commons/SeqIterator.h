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
#include "kseq.h"
#include "KSeqBufferReader.h"
#include "ProdigalWrapper.h"

#ifdef OPENMP
#include <omp.h>
#endif

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define kmerLength 8

#define nuc2int(x) (x & 14u)>>1u

using namespace std;

typedef struct PredictedBlock {
    PredictedBlock(uint32_t start, uint32_t end, int strand) : start(start), end(end), strand(strand) { }
    uint32_t start;
    uint32_t end;
    int strand; //true for forward
}PredictedBlock;

class SeqIterator
{
private:
    string iRCT;
    string atcg;
    vector<int> aaFrames[6];
    uint64_t powers[10];
    int nuc2aa[8][8][8];
    uint64_t nuc2num[4][4][4];

    void addDNAInfo_QueryKmer(uint64_t & kmer, const char * seq, int forOrRev, const int & kmerCnt, const int & frame);
    void addDNAInfo_TargetKmer(uint64_t & kmer, const char * seq, const PredictedBlock& block, const int & startOfKmer);

public:
    void fillQueryKmerBuffer(const char * seq , QueryKmerBuffer & kmerBuffer, size_t & posToWrite, const int & seqID);
    string reverseCompliment(string & read) const ;
    void sixFrameTranslation(const char * seq);
    void translateBlock(const char* seq, PredictedBlock & block);
    void getTranslationBlocks(struct _gene * genes, struct _node * nodes, PredictedBlock * blocks, size_t numOfGene, size_t length, size_t & numOfBlocks);
    void getTranslationBlocks2(struct _gene * genes, struct _node * nodes, vector<PredictedBlock> & blocks, size_t numOfGene, size_t length, size_t & numOfBlocks);


    size_t kmerNumOfSixFrameTranslation(const string & seq);
    size_t getNumOfKmerForBlock(const PredictedBlock & block);
    void fillBufferWithKmerFromBlock(const PredictedBlock & block, const char * seq, TargetKmerBuffer & kmerBuffer, size_t & posToWrite, const int & seqID, int taxIdAtRank);
    void printKmerInDNAsequence(uint64_t kmer);
    SeqIterator();
};
#endif //ADKMER4_KMEREXTRACTOR_H

