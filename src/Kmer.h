//
// Created by KJB on 25/08/2020.
//

#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
typedef struct KmerInfo {
    int sequenceID;
    int pos;
    bool redundancy;
} KmerInfo;

typedef struct Kmer {
    uint64_t ADkmer;
    KmerInfo info;
} Kmer;

typedef struct matchedKmer{
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
