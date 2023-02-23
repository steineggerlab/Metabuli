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
#include "xxh3.h"
#include <queue>
#include "LocalParameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define kmerLength 8

#define nuc2int(x) (x & 14u)>>1u

using namespace std;

typedef struct PredictedBlock {
    PredictedBlock(int start, int end, int strand) : start(start), end(end), strand(strand) {}

    void printPredictedBlock() {
        cout << strand << " " << start << " " << end << endl;
    }

    int start;
    int end;
    int strand; //true for forward
} PredictedBlock;

class SeqIterator {
private:
    static const string iRCT;
    static const string atcg;
    vector<int> aaFrames[6];
    uint64_t powers[10];
    int nuc2aa[8][8][8];
    uint64_t nuc2num[4][4][4];
    uint32_t * mask;
    int * mask_int;
    uint32_t spaceNum;
    int spaceNum_int;
    int bitsForCodon;
    int bitsFor8Codons;

    void addDNAInfo_QueryKmer(uint64_t &kmer, const char *seq, int forOrRev, const int &kmerCnt, const int &frame,
                              int readLength);

    void addDNAInfo_TargetKmer(uint64_t &kmer, const char *seq, const PredictedBlock &block, const int &kmerCnt);

public:
    void fillQueryKmerBuffer(const char *seq, int seqLen, QueryKmerBuffer &kmerBuffer, size_t &posToWrite, const int &seqID,
                             uint32_t offset = 0);

    string reverseCompliment(string &read) const;

    char *reverseCompliment(char *read, int length) const;

    void sixFrameTranslation(const char *seq);

    bool translateBlock(const char *seq, PredictedBlock block);

    void generateIntergenicKmerList(struct _gene *genes, struct _node *nodes, int numberOfGenes,
                                    vector<uint64_t> &intergenicKmerList, const char *seq);

    void getExtendedORFs(struct _gene *genes, struct _node *nodes, vector<PredictedBlock> &blocks, size_t numOfGene,
            size_t length, size_t &numOfBlocks, vector<uint64_t> &intergenicKmerList, const char *seq);

    void getMinHashList(priority_queue<uint64_t> &sortedHashQue, const char *seq);

    bool
    compareMinHashList(priority_queue<uint64_t> list1, priority_queue<uint64_t> &list2, size_t length1, size_t length2);

    static size_t kmerNumOfSixFrameTranslation(const char *seq);

    size_t getNumOfKmerForBlock(const PredictedBlock &block);

    int fillBufferWithKmerFromBlock(const PredictedBlock &block, const char *seq, TargetKmerBuffer &kmerBuffer,
                                     size_t &posToWrite, const uint32_t &seqID, int taxIdAtRank);

    void printKmerInDNAsequence(uint64_t kmer);

    explicit SeqIterator(const LocalParameters &par);
    ~SeqIterator();
};

#endif //ADKMER4_KMEREXTRACTOR_H

