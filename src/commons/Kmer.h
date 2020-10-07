//
// Created by KJB on 25/08/2020.
//

#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
typedef struct KmerInfo {
    KmerInfo(int seqID, int pos, bool redundancy) : sequenceID(seqID), pos(pos), redundancy(redundancy) {}
    int sequenceID;
    int pos;
    bool redundancy = false;
} KmerInfo;

typedef struct Kmer {
    Kmer(uint64_t ADkmer, int seqID, int pos, bool redundacy) : ADkmer(ADkmer),info(seqID, pos, redundacy) {}
    uint64_t ADkmer;
    KmerInfo info;
} Kmer;

typedef struct MatchedKmer{
    MatchedKmer(int quID, int tarID, int taxID, int pos, uint8_t hamming, bool red): queryID(quID), tragetID(tarID), taxID(taxID), posInfo(pos), hammingDistance(hamming), redundancy(red) {}
    int queryID;
    int tragetID;
    int taxID;
    int posInfo;
    uint8_t hammingDistance;
    bool redundancy;
} __attribute__((packed)) matchedKmer;

//typedef struct diffIdxKmer{
//    uint16_t * diffIdx;
//    int sequenceID;
//    int pos;
//} diffIdxKmer;
#endif //ADKMER3_KMER_H
