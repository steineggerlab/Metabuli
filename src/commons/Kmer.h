//
// Created by KJB on 25/08/2020.
//

#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
////
typedef struct KmerInfo {
    KmerInfo(int seqID = 0, uint32_t pos = 0, bool redundancy = 0) : sequenceID(seqID), pos(pos), redundancy(redundancy) {}
    int sequenceID;
    uint32_t pos;
    bool redundancy = false;
} KmerInfo;

typedef struct Kmer {
    Kmer(uint64_t ADkmer, int seqID, uint32_t pos, bool redundacy) : ADkmer(ADkmer),info(seqID, pos, redundacy) {}
    uint64_t ADkmer;
    KmerInfo info;
} Kmer;

typedef struct TargetKmerInfo{
    TargetKmerInfo(int seqID = 0, bool redundancy = 0) : sequenceID(seqID), redundancy(redundancy) {}
    int sequenceID;
    bool redundancy;
} TargetKmerInfo;

typedef struct TargetKmer{
    TargetKmer(uint64_t ADkmer, int seqID, bool redundacy) : ADkmer(ADkmer),info(seqID, redundacy) {}
    uint64_t ADkmer;
    TargetKmerInfo info;
};

typedef struct MatchedKmer{
    MatchedKmer(int quID, int tarID, int taxID, uint32_t pos, uint8_t hamming, bool red): queryID(quID), tragetID(tarID), taxID(taxID), queryPos(pos), hammingDistance(hamming), redundancy(red) {}
    int queryID;
    int tragetID;
    int taxID;
    uint32_t queryPos;
    uint8_t hammingDistance;
    bool redundancy;
} __attribute__((packed)) matchedKmer;

//typedef struct diffIdxKmer{
//    uint16_t * diffIdx;
//    int sequenceID;
//    int pos;
//} diffIdxKmer;
#endif //ADKMER3_KMER_H
