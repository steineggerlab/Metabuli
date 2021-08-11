//
// Created by KJB on 25/08/2020.
//

#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
#include "NcbiTaxonomy.h"

typedef struct QueryKmerInfo {
    QueryKmerInfo(int seqID = 0, uint32_t pos = 0, uint8_t frame = 0 ) : sequenceID(seqID), pos(pos), frame(frame) {}
    uint32_t sequenceID;
    uint16_t pos;
    uint8_t frame; // 0, 1, 2 are forward, and 3, 4, 5 are reverse
} QueryKmerInfo;

typedef struct QueryKmer {
    QueryKmer(uint64_t ADkmer, int seqID, uint32_t pos, uint32_t isReverse) : ADkmer(ADkmer), info(seqID, pos, isReverse) {}
    QueryKmer():ADkmer(0), info(0,0,0){}
    uint64_t ADkmer;
    QueryKmerInfo info;
} QueryKmer;


struct TargetKmerInfo{
    TargetKmerInfo(int seqID = 0, bool redundancy = false) : sequenceID(seqID), redundancy(redundancy) {}
    uint32_t sequenceID;
    bool redundancy;
};

struct TargetKmer{
    TargetKmer():ADkmer(0), taxIdAtRank(0), info(0, false){};
    TargetKmer(uint64_t ADkmer, TaxID taxIdAtRank, uint32_t seqID, bool redundacy) : ADkmer(ADkmer), taxIdAtRank(taxIdAtRank),info(seqID, redundacy) {}
    uint64_t ADkmer;
    TaxID taxIdAtRank;
    TargetKmerInfo info;
};

struct DiffIdxSplit{
    DiffIdxSplit(uint64_t ADkmer, size_t diffIdxOffset, size_t infoIdxOffset) : ADkmer(ADkmer), diffIdxOffset(diffIdxOffset), infoIdxOffset(infoIdxOffset) { }
    DiffIdxSplit(const DiffIdxSplit & copy) {ADkmer = copy.ADkmer; diffIdxOffset = copy.diffIdxOffset; infoIdxOffset=copy.infoIdxOffset;}
    DiffIdxSplit() {};
    uint64_t ADkmer;
    size_t diffIdxOffset;
    size_t infoIdxOffset;
};


#endif //ADKMER3_KMER_H
