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
#include <functional>
//#include "xxh3.h"
#include "xxhash.h"
#include <queue>
#include "LocalParameters.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include <cstdint>

#ifdef OPENMP
#include <omp.h>
#endif

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define kmerLength 8

#define nuc2int(x) (x & 14u)>>1u

using namespace std;

class SeqIterator {
private:
    // vector<int> aaFrames[6];
    uint64_t powers[10];
    int nuc2aa[8][8][8];
    uint64_t nuc2num[4][4][4];
    uint32_t * mask;
    int * mask_int;
    uint32_t spaceNum;
    int spaceNum_int;
    int bitsForCodon;
    int bitsFor8Codons;
    int smerLen;

    void addDNAInfo_QueryKmer(uint64_t &kmer, const char *seq, int forOrRev, uint32_t kmerCnt, uint32_t frame,
                              int readLength);

    // void addDNAInfo_TargetKmer(uint64_t &kmer, const char *seq, const PredictedBlock &block, int kmerCnt);

public:
    static const string iRCT;
    static const string atcg;
    
    string reverseComplement(string &read) const;
    
    void devideToCdsAndNonCds(const char *maskedSeq,
                              size_t seqLen,
                              const vector<CDSinfo> &cdsInfo, 
                              vector<string> &cds,
                              vector<string> &nonCds);
    
    void fillQueryKmerBuffer(const char *seq, int seqLen, QueryKmerBuffer &kmerBuffer, size_t &posToWrite,
                             uint32_t seqID, vector<int> *aaFrames, uint32_t offset = 0);

    void fillQuerySyncmerBuffer(const char *seq, int seqLen, QueryKmerBuffer &kmerBuffer, size_t &posToWrite,
                             uint32_t seqID, vector<int> *aaFrames, uint32_t offset = 0);

    string reverseCompliment(string &read) const;

    char *reverseCompliment(char *read, size_t length) const;

    void sixFrameTranslation(const char *seq, int seqLen, vector<int> *aaFrames);

    bool translateBlock(const char *seq, PredictedBlock block, vector<int> & aaSeq, size_t length);

    void translate(const string & seq, vector<int> & aa, int frame = 0) {
        aa.clear();
        if(aa.capacity() < seq.length() / 3 + 1) {
            aa.reserve(seq.length() / 3 + 1);
        }
        for (int i = 0 + frame; i + 2 < (int) seq.length(); i = i + 3) {
            aa.push_back(nuc2aa[nuc2int(atcg[seq[i]])][nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])]);
        }
    }

    void generateIntergenicKmerList(struct _gene *genes, struct _node *nodes, int numberOfGenes,
                                    vector<uint64_t> &intergenicKmerList, const char *seq);

    void getMinHashList(priority_queue<uint64_t> &sortedHashQue, const char *seq);

    bool
    compareMinHashList(priority_queue<uint64_t> list1, priority_queue<uint64_t> &list2, size_t length1, size_t length2);

    static size_t kmerNumOfSixFrameTranslation(const char *seq);

    size_t getNumOfKmerForBlock(const PredictedBlock &block);

    // int fillBufferWithKmerFromBlock(const PredictedBlock &block,
    //                                 const char *seq, 
    //                                 TargetKmerBuffer &kmerBuffer,
    //                                 size_t &posToWrite, 
    //                                 int seqID, 
    //                                 int taxIdAtRank, 
    //                                 const vector<int> & aaSeq);

    int fillBufferWithKmerFromBlock(const char *seq,
                                    TargetKmerBuffer &kmerBuffer,
                                    size_t &posToWrite,
                                    int seqID,
                                    int taxIdAtRank,
                                    const vector<int> & aaSeq,
                                    int blockStrand = 0,
                                    int blockStart = 0,
                                    int blockEnd = 0);
    
    int fillBufferWithSyncmer(const char *seq,
                              TargetKmerBuffer &kmerBuffer,
                              size_t &posToWrite,
                              int seqID,
                              int taxIdAtRank,
                              const vector<int> & aaSeq,
                              int blockStrand = 0,
                              int blockStart = 0,
                              int blockEnd = 0);

    bool isSyncmer(const vector<int> &aaSeq, int startPos, int k, int s) {
        int min_smer_value = INT32_MAX;
        int min_smer_pos = -1;

        for (int i = 0; i <= k - s; ++i) {
            int current_value = 0;
            for (int j = 0; j < s; ++j) {
                current_value += aaSeq[startPos + i + j] * powers[j];
            }
            if (current_value < min_smer_value) {
                min_smer_value = current_value;
                min_smer_pos = i;
            }
        }
        return (min_smer_pos == 0 || min_smer_pos == (k - s));
    }
    // bool isSyncmer(const vector<int> & aaSeq, int startPos, int k, int s);

    void addDNAInfo_TargetKmer(uint64_t & kmer, const char * seq, int kmerCnt, int strand, int start, int end);

    static void maskLowComplexityRegions(const char * seq, char * maskedSeq, ProbabilityMatrix & probMat,
                                         float maskProb, const BaseMatrix * subMat);

    void printKmerInDNAsequence(uint64_t kmer);

    void printAAKmer(uint64_t kmer, int shits = 28);

    explicit SeqIterator(const LocalParameters &par);
    ~SeqIterator();
};

#endif //ADKMER4_KMEREXTRACTOR_H

