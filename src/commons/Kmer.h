//
// Created by KJB on 25/08/2020.
//

#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
////
typedef struct QueryKmerInfo {
    QueryKmerInfo(int seqID = 0, uint32_t pos = 0, uint32_t isReverse = 0 ) : sequenceID(seqID), pos(pos), frame(isReverse) {}
    int sequenceID;
    uint32_t pos;
    uint32_t frame; // 0, 1, 2 are forward, and 3, 4, 5 are reverse
} QueryKmerInfo;

typedef struct QueryKmer {
    QueryKmer(uint64_t ADkmer, int seqID, uint32_t pos, uint32_t isReverse) : ADkmer(ADkmer), info(seqID, pos, isReverse) {}
    uint64_t ADkmer;
    QueryKmerInfo info;
} QueryKmer;


typedef struct TargetKmerInfo{
    TargetKmerInfo(int seqID = 0, bool redundancy = false) : sequenceID(seqID), redundancy(redundancy) {}
    int sequenceID;
    bool redundancy;
} TargetKmerInfo;

typedef struct TargetKmer{
    TargetKmer(uint64_t ADkmer, int seqID, bool redundacy) : ADkmer(ADkmer),info(seqID, redundacy) {}
    uint64_t ADkmer;
    TargetKmerInfo info;
} TargetKmer;

typedef struct MatchedKmer{
    MatchedKmer(int quID, int tarID, int taxID, uint32_t pos, uint8_t hamming, bool red, int queryStrand): queryID(quID), targetID(tarID), taxID(taxID), queryPos(pos), hammingDistance(hamming), redundancy(red), queryFrame(queryStrand) {}
    int queryID;
    int targetID;
    int taxID;
    uint32_t queryPos;
    uint8_t hammingDistance;
    bool redundancy;
    int queryFrame;
} __attribute__((packed)) matchedKmer;

//typedef struct diffIdxKmer{
//    uint16_t * diffIdx;
//    int sequenceID;
//    int pos;
//} diffIdxKmer;
#endif //ADKMER3_KMER_H
