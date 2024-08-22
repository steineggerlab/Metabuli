#ifndef METABULI_KMEREXTRACTER_H
#define METABULI_KMEREXTRACTER_H
#include "SeqIterator.h"
#include "QueryIndexer.h"
#include "KSeqWrapper.h"
#include "common.h"
#include <unordered_map>

class KmerExtractor {
private:
    SeqIterator * seqIterator;
    // Parameters
    int spaceNum;
    int maskMode;
    float maskProb;

    // For masking reads
    ProbabilityMatrix * probMatrix;
    BaseMatrix * subMat;

    // Extract query k-mer
    void fillQueryKmerBufferParallel(KSeqWrapper* kseq1,
                                     QueryKmerBuffer &kmerBuffer,
                                     vector<Query> & queryList,
                                     const QuerySplit & currentSplit,
                                     const LocalParameters &par);

    void fillQueryKmerBufferParallel_paired(KSeqWrapper* kseq1,
                                            KSeqWrapper* kseq2,
                                            QueryKmerBuffer &kmerBuffer,
                                            vector<Query> &queryList,
                                            const QuerySplit & currentSplit,
                                            const LocalParameters &par);

    void fillQueryKmerBufferParallel_paired2(KSeqWrapper* kseq1,
                                            KSeqWrapper* kseq2,
                                            QueryKmerBuffer &kmerBuffer,
                                            vector<Query> &queryList,
                                            const QuerySplit & currentSplit,
                                            const LocalParameters &par);                                        

public:
    explicit KmerExtractor(const LocalParameters & par);
    ~KmerExtractor();
    void extractQueryKmers(QueryKmerBuffer &kmerBuffer,
                           vector<Query> & queryList,
                           const QuerySplit & currentSplit,
                           const LocalParameters &par,
                           KSeqWrapper* kseq1,
                           KSeqWrapper* kseq2 = nullptr);


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
