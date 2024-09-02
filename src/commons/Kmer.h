#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
#include "NcbiTaxonomy.h"

struct QueryKmerInfo {
    explicit QueryKmerInfo(uint32_t seqID = 0, uint32_t pos = 0, uint8_t frame = 0 ) : pos(pos), sequenceID(seqID), frame(frame) {}
    uint64_t pos : 32;
    uint64_t sequenceID : 29;
    uint64_t frame : 3; // 0, 1, 2 are forward, and 3, 4, 5 are reverse 1 byte
}; // 8 byte


typedef struct QueryKmer {
    QueryKmer(uint64_t ADkmer, uint32_t seqID, uint32_t pos, uint8_t frame) : ADkmer(ADkmer), info(seqID, pos, frame) {}
    QueryKmer():ADkmer(0), info(0,0,0){}
    uint64_t ADkmer; // 8 byte
    QueryKmerInfo info; // 8 byte
} QueryKmer; // 16 byte


struct TargetKmerInfo{
    explicit TargetKmerInfo(int seqID = 0, bool redundancy = false) : sequenceID(seqID), redundancy(redundancy) {}
    int sequenceID : 31;
    int redundancy : 1;
    bool operator == (const TargetKmerInfo & info) const{
        return (sequenceID == info.sequenceID && this->redundancy==info.redundancy);
    }
}; // 4 bytes

struct TargetKmer{
    TargetKmer(): info(0, false), taxIdAtRank(0), ADkmer(0)  { };
    TargetKmer(uint64_t ADkmer, TaxID taxIdAtRank, int seqID, bool redundacy)
        : info(seqID, redundacy), taxIdAtRank(taxIdAtRank), ADkmer(ADkmer) {}
    TargetKmerInfo info; // 4 byte
    TaxID taxIdAtRank; // 4 byte
    uint64_t ADkmer; // 8 byte
};

struct DiffIdxSplit{
    DiffIdxSplit(uint64_t ADkmer, size_t diffIdxOffset, size_t infoIdxOffset) : ADkmer(ADkmer), diffIdxOffset(diffIdxOffset), infoIdxOffset(infoIdxOffset) { }
    DiffIdxSplit(const DiffIdxSplit & copy) {ADkmer = copy.ADkmer; diffIdxOffset = copy.diffIdxOffset; infoIdxOffset=copy.infoIdxOffset;}
    DiffIdxSplit() {};
    DiffIdxSplit& operator=(const DiffIdxSplit&) = default;
    uint64_t ADkmer;
    size_t diffIdxOffset;
    size_t infoIdxOffset;
};


#endif //ADKMER3_KMER_H
