#ifndef METABULI_KMEREXTRACTER_H
#define METABULI_KMEREXTRACTER_H
#include "SeqIterator.h"
#include "QueryIndexer.h"
#include "KSeqWrapper.h"
#include "common.h"
#include <unordered_map>
#include <atomic>
#include <cstdint>

class KmerExtractor {
private:
    const LocalParameters &par;
    KmerScanner ** kmerScanners;
    
    // Parameters
    int spaceNum;
    int maskMode;
    float maskProb;

    // For masking reads
    ProbabilityMatrix * probMatrix;
    BaseMatrix * subMat;

    // Extract query k-mer
    void fillQueryKmerBufferParallel(KSeqWrapper* kseq1,
                                     Buffer<QueryKmer> &kmerBuffer,
                                     vector<Query> & queryList,
                                     const QuerySplit & currentSplit,
                                     const LocalParameters &par);

    void fillQueryKmerBufferParallel_paired(KSeqWrapper* kseq1,
                                            KSeqWrapper* kseq2,
                                            Buffer<QueryKmer> &kmerBuffer,
                                            vector<Query> &queryList,
                                            const QuerySplit & currentSplit,
                                            const LocalParameters &par);

    void loadChunkOfReads(KSeqWrapper *kseq,
                          vector<Query> & queryList,
                          size_t & processedQueryNum,
                          size_t chunkSize,
                          size_t chunkEnd,
                          vector<string> & reads,
                          vector<bool> & emptyReads,
                          size_t & count,
                          bool isReverse);

    void processSequence(
        size_t count,
        size_t processedQueryNum,
        const vector<string> & reads,
        const vector<bool> & emptyReads,
        char *seq,
        char *maskedSeq,
        size_t & maxReadLength,
        Buffer<QueryKmer> &kmerBuffer,
        const vector<Query> & queryList,
        bool isReverse);

    void fillQueryKmerBuffer(
        const char *seq,
        int seqLen, 
        Buffer<QueryKmer> &kmerBuffer, 
        size_t &posToWrite, 
        uint32_t seqID, 
        uint32_t offset = 0);
                                      
public:
    explicit KmerExtractor(
        const LocalParameters & par,
        const GeneticCode &geneticCode,
        int kmerFormat);
    ~KmerExtractor();
    void extractQueryKmers(Buffer<QueryKmer> &kmerBuffer,
                           vector<Query> & queryList,
                           const QuerySplit & currentSplit,
                           const LocalParameters &par,
                           KSeqWrapper* kseq1,
                           KSeqWrapper* kseq2 = nullptr);

    int extractTargetKmers(
        const char *seq,
        Buffer<Kmer_union> &kmerBuffer,
        size_t &posToWrite,
        int seqID,
        int taxIdAtRank,
        SequenceBlock block);
};

static inline bool compareForLinearSearch(const QueryKmer &a, const QueryKmer &b) {
    if (a.ADkmer < b.ADkmer) {
        return true;
    } else if (a.ADkmer == b.ADkmer) {
        return (a.info.sequenceID < b.info.sequenceID);
    }
    return false;
}

#endif //METABULI_KMEREXTRACTER_H
